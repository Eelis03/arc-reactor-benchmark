"""Zero-dimensional plasma power balance.

The plasma is described by one volume-averaged density, one volume-averaged
temperature shared by electrons and ions, a composition, and a shape. Every term
in the balance is a function of those numbers and of the geometry.

The balance solved elsewhere in this package is

    P_alpha + P_aux = P_radiated + P_transport

where the transport term is the stored thermal energy divided by the confinement
time. This module computes each term. It does not solve anything.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Final

from arc_benchmark.model.confinement import ConfinementInputs
from arc_benchmark.model.constants import (
    ALPHA_FRACTION,
    DT_ENERGY_KEV,
    JOULE_PER_KEV,
    NEUTRON_FRACTION,
)
from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.model.profiles import FLAT_PROFILE, ProfileShape
from arc_benchmark.model.radiation import (
    ImpurityRadiator,
    bremsstrahlung_density,
    line_radiation_density,
    synchrotron_density,
)
from arc_benchmark.model.reactivity import dt_reactivity

__all__ = [
    "LossPowerConvention",
    "PlasmaComposition",
    "PlasmaState",
    "PowerTerms",
    "power_terms",
]

_WATT_PER_MEGAWATT: Final[float] = 1.0e6
_JOULE_PER_MEGAJOULE: Final[float] = 1.0e6


class LossPowerConvention(Enum):
    """Which power the confinement scaling is evaluated at.

    The choice matters and there is no universally correct answer, so it is a
    parameter rather than a hidden decision.

    ``SEPARATRIX`` treats core radiation as a loss channel parallel to transport
    and evaluates the scaling at the power actually crossing the last closed flux
    surface, ``P_alpha + P_aux - P_rad``. This is what reactor systems studies
    normally do, and it is the only convention under which the radiation terms
    appear explicitly in the balance.

    ``TOTAL`` evaluates the scaling at ``P_alpha + P_aux``, which is the
    definition the IPB98 regression itself was fitted with: the ITER Physics
    Basis defines the loss power as total input power less the rate of change of
    stored energy, without subtracting radiation. Under this convention the
    radiated power is reported but is already inside the transport channel and is
    not subtracted again.

    The two differ by the radiated fraction, which for an ARC-class point in this
    model is around a fifth of the heating power, so the choice is worth a
    visible amount and is swept in the results.
    """

    SEPARATRIX = "separatrix"
    TOTAL = "total"


@dataclass(frozen=True, slots=True)
class PlasmaComposition:
    """Fuel dilution, effective charge, and the impurities that cause both.

    Quasineutrality and the definition of the effective charge are enforced
    rather than assumed, so a composition cannot report a fuel fraction that is
    inconsistent with its impurity content.

    Attributes:
        helium_fraction: Thermalised helium ash density over electron density.
        impurities: Impurity species present, each with its charge, its
            concentration, and its radiative loss parameter.
        mass_number: Average hydrogenic isotope mass in atomic mass units.
    """

    helium_fraction: float = 0.0
    impurities: tuple[ImpurityRadiator, ...] = ()
    mass_number: float = 2.5

    def __post_init__(self) -> None:
        """Reject compositions that leave no fuel or that are not physical."""
        if not 0.0 <= self.helium_fraction < 0.5:
            raise ValueError(f"helium_fraction must lie in [0, 0.5), got {self.helium_fraction}")
        if self.mass_number <= 0.0:
            raise ValueError(f"mass_number must be positive, got {self.mass_number}")
        if self.fuel_fraction <= 0.0:
            raise ValueError(
                "impurity and helium content leaves no fuel: fuel fraction is "
                f"{self.fuel_fraction}"
            )

    @property
    def fuel_fraction(self) -> float:
        """``n_DT / n_e``, from quasineutrality with fully stripped impurities."""
        impurity_charge = sum(imp.atomic_number * imp.concentration for imp in self.impurities)
        return 1.0 - 2.0 * self.helium_fraction - impurity_charge

    @property
    def z_effective(self) -> float:
        """``sum(n_j Z_j**2) / n_e``, the effective charge.

        Reduces to ``1 + 2 f_He + sum(c_z Z (Z - 1))``, so a pure hydrogenic
        plasma returns exactly one.
        """
        impurity_term = sum(
            imp.concentration * imp.atomic_number * (imp.atomic_number - 1)
            for imp in self.impurities
        )
        return 1.0 + 2.0 * self.helium_fraction + impurity_term

    @property
    def ion_fraction(self) -> float:
        """``n_i / n_e``, the total ion density over electron density.

        Below one for any plasma containing an impurity, which is why the
        pressure of a diluted plasma is below twice ``n_e T``.
        """
        return (
            self.fuel_fraction
            + self.helium_fraction
            + sum(imp.concentration for imp in self.impurities)
        )


@dataclass(frozen=True, slots=True)
class PlasmaState:
    """One candidate operating point, before any balance has been solved.

    Attributes:
        geometry: Plasma shape.
        toroidal_field: Vacuum toroidal field on axis in tesla.
        plasma_current_ma: Plasma current in megaampere.
        electron_density: Volume-averaged electron density in inverse cubic metre.
        temperature_kev: Volume-averaged temperature in keV, shared by electrons
            and ions.
        composition: Fuel dilution and impurity content.
        profile: Profile shape used for the volume-integral correction factors.
            The default is flat, which makes every factor exactly one and the
            model genuinely zero-dimensional.
        confinement_multiplier: The H factor multiplying the scaling law.
        wall_reflectivity: Fraction of synchrotron emission reflected back into
            the plasma by the first wall.
    """

    geometry: PlasmaGeometry
    toroidal_field: float
    plasma_current_ma: float
    electron_density: float
    temperature_kev: float
    composition: PlasmaComposition = field(default_factory=PlasmaComposition)
    profile: ProfileShape = FLAT_PROFILE
    confinement_multiplier: float = 1.0
    wall_reflectivity: float = 0.9

    def __post_init__(self) -> None:
        """Reject states no balance could be evaluated at."""
        for name in ("toroidal_field", "plasma_current_ma", "electron_density", "temperature_kev"):
            value = getattr(self, name)
            if value <= 0.0:
                raise ValueError(f"{name} must be positive, got {value}")
        if self.confinement_multiplier <= 0.0:
            raise ValueError(
                f"confinement_multiplier must be positive, got {self.confinement_multiplier}"
            )

    @property
    def ion_density(self) -> float:
        """Total ion density in inverse cubic metre."""
        return self.composition.ion_fraction * self.electron_density

    @property
    def fuel_density(self) -> float:
        """Deuterium plus tritium density in inverse cubic metre."""
        return self.composition.fuel_fraction * self.electron_density

    @property
    def thermal_pressure_pa(self) -> float:
        """Volume-averaged thermal pressure in pascal.

        ``<p> = (n_e + n_i) T e_keV`` with the profile correction for
        ``<n T> / (<n> <T>)`` applied. Fast alpha pressure is not included; a
        design point where it matters is outside what this model supports.
        """
        density_sum = self.electron_density + self.ion_density
        return (
            density_sum
            * self.temperature_kev
            * JOULE_PER_KEV
            * self.profile.stored_energy_factor()
        )

    @property
    def stored_energy_mj(self) -> float:
        """Thermal stored energy in megajoule, ``(3/2) <p> V``."""
        return 1.5 * self.thermal_pressure_pa * self.geometry.volume / _JOULE_PER_MEGAJOULE

    def confinement_inputs(self, loss_power_mw: float) -> ConfinementInputs:
        """Assemble the inputs a confinement scaling needs, at a given loss power.

        The line-averaged density the scalings were fitted against is taken equal
        to the volume-averaged density. For a peaked profile the line average is
        higher, so this is a conservative choice, and it is stated in the design
        notes rather than corrected with a factor this model cannot justify.
        """
        return ConfinementInputs(
            plasma_current_ma=self.plasma_current_ma,
            toroidal_field=self.toroidal_field,
            line_averaged_density_e19=self.electron_density / 1.0e19,
            loss_power_mw=loss_power_mw,
            major_radius=self.geometry.major_radius,
            inverse_aspect_ratio=self.geometry.inverse_aspect_ratio,
            areal_elongation=self.geometry.areal_elongation,
            mass_number=self.composition.mass_number,
        )


@dataclass(frozen=True, slots=True)
class PowerTerms:
    """Every source and sink in the balance, in megawatt, at a fixed state.

    The auxiliary power is not here: it is the unknown the operating point solver
    determines. Everything in this object is a function of the plasma state
    alone.

    Attributes:
        fusion_power_mw: Total deuterium tritium fusion power.
        alpha_power_mw: Fusion power deposited in the plasma.
        neutron_power_mw: Fusion power leaving with the neutrons.
        bremsstrahlung_mw: Bremsstrahlung loss.
        synchrotron_mw: Synchrotron loss after wall reflection.
        line_radiation_mw: Impurity line radiation loss.
        stored_energy_mj: Thermal stored energy in megajoule.
        reactivity_m3_s: The Maxwellian reactivity used, in cubic metre per second.
    """

    fusion_power_mw: float
    alpha_power_mw: float
    neutron_power_mw: float
    bremsstrahlung_mw: float
    synchrotron_mw: float
    line_radiation_mw: float
    stored_energy_mj: float
    reactivity_m3_s: float

    @property
    def radiated_power_mw(self) -> float:
        """Total radiated power in megawatt."""
        return self.bremsstrahlung_mw + self.synchrotron_mw + self.line_radiation_mw


def power_terms(state: PlasmaState) -> PowerTerms:
    """Evaluate every state-dependent term of the power balance.

    Fusion power uses ``n_D n_T <sigma v> E_fus V`` with a fifty-fifty fuel
    mixture, so ``n_D = n_T = n_DT / 2`` and the density product is
    ``(n_DT / 2)**2``. Writing it that way rather than as ``n_DT**2 / 4`` is the
    same arithmetic, but it makes the factor of four visible, and forgetting it is
    the single most common error in a zero-dimensional fusion power estimate.

    Args:
        state: The plasma state to evaluate.

    Returns:
        Every source and sink except auxiliary heating.
    """
    volume = state.geometry.volume
    reactivity = float(dt_reactivity(state.temperature_kev))

    fuel_half = 0.5 * state.fuel_density
    fusion_density = (
        fuel_half**2 * reactivity * DT_ENERGY_KEV * JOULE_PER_KEV
    ) * state.profile.fusion_factor(state.temperature_kev)
    fusion_mw = fusion_density * volume / _WATT_PER_MEGAWATT

    bremsstrahlung_mw = (
        float(
            bremsstrahlung_density(
                state.electron_density,
                state.temperature_kev,
                state.composition.z_effective,
            )
        )
        * state.profile.bremsstrahlung_factor()
        * volume
        / _WATT_PER_MEGAWATT
    )
    synchrotron_mw = (
        float(
            synchrotron_density(
                state.electron_density,
                state.temperature_kev,
                state.toroidal_field,
                state.geometry.minor_radius,
                state.wall_reflectivity,
            )
        )
        * state.profile.synchrotron_factor()
        * volume
        / _WATT_PER_MEGAWATT
    )
    line_mw = (
        float(line_radiation_density(state.electron_density, state.composition.impurities))
        * state.profile.density_square_factor()
        * volume
        / _WATT_PER_MEGAWATT
    )

    return PowerTerms(
        fusion_power_mw=fusion_mw,
        alpha_power_mw=ALPHA_FRACTION * fusion_mw,
        neutron_power_mw=NEUTRON_FRACTION * fusion_mw,
        bremsstrahlung_mw=bremsstrahlung_mw,
        synchrotron_mw=synchrotron_mw,
        line_radiation_mw=line_mw,
        stored_energy_mj=state.stored_energy_mj,
        reactivity_m3_s=reactivity,
    )
