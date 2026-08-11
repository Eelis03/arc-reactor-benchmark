"""Energy confinement time scaling laws.

Three empirical scalings are provided. They are written in one canonical basis,
``(I_p, B_t, n_e19, P_loss, R_0, epsilon, kappa_a, M)``, so that the same solver
and the same tests apply to all of them, and the published form of each is
recorded in its docstring so the translation into that basis can be checked.

None of these functions knows anything about a power balance. They map plasma
parameters to a confinement time and nothing else, which is what allows the
algorithm layer to accept any of them through a structural Protocol without
importing this module.

A scaling law extrapolated a factor of two outside the database it was fitted to
is the dominant uncertainty in any projection of this kind. That is why more than
one is provided: every conclusion this package draws is reported against all
three so that its dependence on the choice is visible rather than hidden.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

__all__ = [
    "CONFINEMENT_SCALINGS",
    "IPB98Y2",
    "ITER89P",
    "PETTY08",
    "ConfinementInputs",
    "PowerLawConfinement",
    "ScalingExponents",
]


@dataclass(frozen=True, slots=True)
class ConfinementInputs:
    """The plasma parameters every scaling in this module is a function of.

    Attributes:
        plasma_current_ma: Plasma current in megaampere.
        toroidal_field: Vacuum toroidal field on axis in tesla.
        line_averaged_density_e19: Line-averaged electron density in units of
            1e19 inverse cubic metre.
        loss_power_mw: Power crossing the separatrix in megawatt.
        major_radius: Major radius in metre.
        inverse_aspect_ratio: ``a / R``, dimensionless.
        areal_elongation: ``kappa_a = S / (pi a**2)``, dimensionless.
        mass_number: Average hydrogenic isotope mass in atomic mass units. 2.5
            for a fifty-fifty deuterium tritium mixture.
    """

    plasma_current_ma: float
    toroidal_field: float
    line_averaged_density_e19: float
    loss_power_mw: float
    major_radius: float
    inverse_aspect_ratio: float
    areal_elongation: float
    mass_number: float = 2.5

    def __post_init__(self) -> None:
        """Reject inputs a power law cannot be evaluated at."""
        for field_name in (
            "plasma_current_ma",
            "toroidal_field",
            "line_averaged_density_e19",
            "loss_power_mw",
            "major_radius",
            "areal_elongation",
            "mass_number",
        ):
            value = getattr(self, field_name)
            if value <= 0.0:
                raise ValueError(f"{field_name} must be positive, got {value}")
        if not 0.0 < self.inverse_aspect_ratio < 1.0:
            raise ValueError(
                f"inverse_aspect_ratio must lie in (0, 1), got {self.inverse_aspect_ratio}"
            )

    def with_loss_power(self, loss_power_mw: float) -> ConfinementInputs:
        """Return a copy at a different loss power, leaving everything else alone."""
        return ConfinementInputs(
            plasma_current_ma=self.plasma_current_ma,
            toroidal_field=self.toroidal_field,
            line_averaged_density_e19=self.line_averaged_density_e19,
            loss_power_mw=loss_power_mw,
            major_radius=self.major_radius,
            inverse_aspect_ratio=self.inverse_aspect_ratio,
            areal_elongation=self.areal_elongation,
            mass_number=self.mass_number,
        )


@dataclass(frozen=True, slots=True)
class ScalingExponents:
    """Exponents of a confinement scaling in the canonical basis."""

    current: float
    field: float
    density: float
    power: float
    major_radius: float
    inverse_aspect_ratio: float
    elongation: float
    mass: float


@dataclass(frozen=True, slots=True)
class PowerLawConfinement:
    """A confinement scaling that is a product of powers of the canonical inputs.

    Attributes:
        name: Short identifier used in reports and traces.
        coefficient: Multiplicative constant in the canonical basis.
        exponents: Exponents in the canonical basis.
        confinement_mode: ``"H"`` or ``"L"``, the regime the fit was taken in.
        published_form: The scaling exactly as printed in its source.
        reference: Author, venue, and year of the source.
    """

    name: str
    coefficient: float
    exponents: ScalingExponents
    confinement_mode: str
    published_form: str
    reference: str

    @property
    def power_degradation(self) -> float:
        """Positive ``alpha_P`` such that ``tau_E`` goes as ``P**(-alpha_P)``.

        The operating point solver inverts the balance analytically using this
        number, so it is part of the interface rather than an implementation
        detail.
        """
        return -self.exponents.power

    def tau_e(self, inputs: ConfinementInputs) -> float:
        """Energy confinement time in second at the given plasma parameters."""
        # float(...) because a float raised to a float power is typed as Any: the
        # general case can produce a complex result, which this one never does.
        return float(self.tau_e_at_unit_power(inputs) * inputs.loss_power_mw**self.exponents.power)

    def tau_e_at_unit_power(self, inputs: ConfinementInputs) -> float:
        """Confinement time in second at a loss power of one megawatt.

        Factoring the power dependence out is what makes the steady-state balance
        solvable in closed form: the balance becomes a single power law in the
        loss power, which inverts without iteration.
        """
        e = self.exponents
        return float(
            self.coefficient
            * inputs.plasma_current_ma**e.current
            * inputs.toroidal_field**e.field
            * inputs.line_averaged_density_e19**e.density
            * inputs.major_radius**e.major_radius
            * inputs.inverse_aspect_ratio**e.inverse_aspect_ratio
            * inputs.areal_elongation**e.elongation
            * inputs.mass_number**e.mass
        )


IPB98Y2: Final[PowerLawConfinement] = PowerLawConfinement(
    name="IPB98(y,2)",
    coefficient=0.0562,
    exponents=ScalingExponents(
        current=0.93,
        field=0.15,
        density=0.41,
        power=-0.69,
        major_radius=1.97,
        inverse_aspect_ratio=0.58,
        elongation=0.78,
        mass=0.19,
    ),
    confinement_mode="H",
    published_form=(
        "tau_E,th = 0.0562 I_p^0.93 B_t^0.15 n_e19^0.41 P_loss^-0.69 "
        "R^1.97 eps^0.58 kappa_a^0.78 M^0.19"
    ),
    reference=(
        "ITER Physics Basis, Chapter 2, Nuclear Fusion 39 (1999) 2175, "
        "equation (20), the ELMy H-mode thermal energy confinement scaling."
    ),
)
"""The ELMy H-mode scaling, in the units the ITER Physics Basis publishes it in.

Current in megaampere, field in tesla, line-averaged density in 1e19 inverse
cubic metre, loss power in megawatt, major radius in metre, and confinement time
in second. No translation is needed: the canonical basis of this module is the
published basis for this scaling.
"""

ITER89P: Final[PowerLawConfinement] = PowerLawConfinement(
    name="ITER89-P",
    # The published form carries density in 1e20 inverse cubic metre. Converting
    # to the 1e19 basis multiplies the coefficient by 10**-0.1.
    coefficient=0.048 * 10.0**-0.1,
    exponents=ScalingExponents(
        current=0.85,
        field=0.20,
        density=0.10,
        power=-0.50,
        # Published as R^1.2 a^0.3. With a = eps R this is R^1.5 eps^0.3.
        major_radius=1.5,
        inverse_aspect_ratio=0.30,
        elongation=0.50,
        mass=0.50,
    ),
    confinement_mode="L",
    published_form=(
        "tau_E = 0.048 I_p^0.85 R^1.2 a^0.3 kappa^0.5 n_e20^0.1 B_t^0.2 M^0.5 P_loss^-0.5"
    ),
    reference=(
        "P. N. Yushmanov et al., 'Scalings for tokamak energy confinement', "
        "Nuclear Fusion 30 (1990) 1999. The ITER89-P L-mode power law."
    ),
)
"""L-mode scaling, kept as the pessimistic bound on any projection.

An ARC-class point evaluated in L-mode rather than H-mode is the honest answer to
the question of what happens if the H-mode assumption does not hold, and it is
included so that question can be answered numerically rather than asserted.
"""

PETTY08: Final[PowerLawConfinement] = PowerLawConfinement(
    name="Petty08",
    coefficient=0.052,
    exponents=ScalingExponents(
        current=0.75,
        field=0.30,
        density=0.32,
        power=-0.47,
        major_radius=2.09,
        inverse_aspect_ratio=0.84,
        elongation=0.88,
        mass=0.0,
    ),
    confinement_mode="H",
    published_form=(
        "tau_E = 0.052 I_p^0.75 B_t^0.30 n_e19^0.32 P_loss^-0.47 R^2.09 eps^0.84 kappa_a^0.88"
    ),
    reference=(
        "C. C. Petty, 'Sizing up plasmas using dimensionless parameters', "
        "Physics of Plasmas 15 (2008) 080501, the constant-beta constrained "
        "H-mode fit."
    ),
)
"""H-mode scaling constrained to be dimensionally consistent at constant beta.

The interesting difference from IPB98(y,2) is the weaker current exponent, 0.75
against 0.93, and the weaker power degradation, 0.47 against 0.69. Both push the
projection for a high-current compact device in the same direction, so this
scaling is the one that most changes the answer for an ARC-class point.
"""

CONFINEMENT_SCALINGS: Final[dict[str, PowerLawConfinement]] = {
    scaling.name: scaling for scaling in (IPB98Y2, ITER89P, PETTY08)
}
"""Every scaling in this module, keyed by name, for sweeps over the choice."""


def _check_module_consistency() -> None:
    """Guard against a mistyped exponent silently changing every result.

    Called once at import. The three power degradation exponents are the numbers
    the operating point solver divides by, so a value at or above one would make
    the steady-state balance have no solution rather than raise later inside a
    sweep.
    """
    for scaling in CONFINEMENT_SCALINGS.values():
        if not 0.0 < scaling.power_degradation < 1.0:
            raise ValueError(
                f"{scaling.name} has power degradation {scaling.power_degradation}, which must "
                "lie in (0, 1) for a steady-state balance to exist"
            )
        if not math.isfinite(scaling.coefficient) or scaling.coefficient <= 0.0:
            raise ValueError(f"{scaling.name} has a non-positive coefficient")


_check_module_consistency()
