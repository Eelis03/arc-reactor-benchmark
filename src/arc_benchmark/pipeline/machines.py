"""Published design points, transcribed with their sources.

Three machines are defined. ARC is the design this package is named after, a
compact high-field pilot plant with demountable rare earth barium copper oxide
magnets. ITER is the conventional low-field reference that the confinement
scaling used here was extrapolated for. SPARC is the intermediate step between
them, from the same group as ARC, and it is included because it is the point
where the high-field argument gets tested experimentally.

ARC here is the Affordable, Robust, Compact tokamak of Sorbom and colleagues
(2015). It is a peer-reviewed engineering design for a fusion pilot plant. It is
not the fictional device of a similar nickname, and nothing in this package has
any relationship to that.

Every published number carried in a :class:`PublishedValues` is a design target
from the cited paper, not a measurement. ITER and SPARC have not operated at
these points and ARC has not been built. The benchmark compares this model
against those published projections, which is a check on the model and not a
check on the machines.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Final

from arc_benchmark.algorithm.balance import PlasmaComposition, PlasmaState
from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.model.plant import PlantParameters
from arc_benchmark.model.radiation import CORONAL_LOSS_PARAMETER, ImpurityRadiator

__all__ = ["MACHINES", "MachineCase", "PublishedValues", "machine"]


@dataclass(frozen=True, slots=True)
class PublishedValues:
    """Design values quoted in the source publication.

    Every field is optional because the sources do not all quote the same
    quantities, and inventing a number to fill a gap would defeat the purpose of
    the comparison.

    Attributes:
        fusion_power_mw: Total fusion power.
        fusion_gain: Plasma gain ``Q``.
        confinement_time_s: Energy confinement time.
        normalised_beta: ``beta_N`` in percent metre tesla per megaampere.
        bootstrap_fraction: Fraction of the plasma current self-driven.
        greenwald_fraction: ``n_e / n_G``.
        auxiliary_power_mw: Heating and current drive power coupled to the plasma.
        thermal_power_mw: Heat delivered to the power conversion system.
        net_electric_mw: Net electrical output.
    """

    fusion_power_mw: float | None = None
    fusion_gain: float | None = None
    confinement_time_s: float | None = None
    normalised_beta: float | None = None
    bootstrap_fraction: float | None = None
    greenwald_fraction: float | None = None
    auxiliary_power_mw: float | None = None
    thermal_power_mw: float | None = None
    net_electric_mw: float | None = None


@dataclass(frozen=True, slots=True)
class MachineCase:
    """One published design point, ready to be solved.

    Attributes:
        name: Short name used in reports.
        description: One sentence on what the machine is.
        source: Citation for every published value carried here.
        state: The plasma state as this model represents the design point.
        published: The design values to compare against.
        plant: Balance of plant parameters, or ``None`` for a device that has no
            power conversion system.
        notes: Statements about known differences between the model and the
            source, carried into the benchmark result.
    """

    name: str
    description: str
    source: str
    state: PlasmaState
    published: PublishedValues
    plant: PlantParameters | None = None
    notes: tuple[str, ...] = field(default_factory=tuple)


_ARC = MachineCase(
    name="ARC",
    description=(
        "Affordable, Robust, Compact: a compact high-field fusion pilot plant with "
        "demountable rare earth barium copper oxide toroidal field magnets and a molten "
        "salt liquid immersion blanket."
    ),
    source=(
        "B. N. Sorbom et al., 'ARC: A compact, high-field, fusion nuclear science "
        "facility and demonstration power plant with demountable magnets', Fusion "
        "Engineering and Design 100 (2015) 378."
    ),
    state=PlasmaState(
        geometry=PlasmaGeometry(
            major_radius=3.3, minor_radius=1.13, elongation=1.84, triangularity=0.375
        ),
        toroidal_field=9.2,
        plasma_current_ma=7.8,
        electron_density=1.3e20,
        temperature_kev=14.0,
        composition=PlasmaComposition(helium_fraction=0.05, mass_number=2.5),
        confinement_multiplier=1.8,
        wall_reflectivity=0.9,
    ),
    published=PublishedValues(
        fusion_power_mw=525.0,
        fusion_gain=13.6,
        normalised_beta=2.59,
        bootstrap_fraction=0.63,
        greenwald_fraction=0.67,
        auxiliary_power_mw=38.6,
        thermal_power_mw=708.0,
        net_electric_mw=190.0,
    ),
    plant=PlantParameters(
        thermal_efficiency=0.40,
        blanket_multiplication=1.30,
        heating_wallplug_efficiency=0.55,
        cryoplant_mw=10.0,
        balance_of_plant_fraction=0.04,
        tritium_and_auxiliary_mw=5.0,
    ),
    notes=(
        "The design point is quoted at an H factor of 1.8 relative to IPB98(y,2), well "
        "above the value the scaling was fitted at, and the confinement time therefore "
        "rests on an assumption rather than on the scaling.",
        "The source solves transport with profiles. This model is volume-averaged, so "
        "the flat-profile fusion power is expected to fall short of the published value.",
        "The blanket energy multiplication is not resolved here. 1.30 is taken as "
        "representative of a lithium blanket with a tritium breeding ratio near 1.1, "
        "where the exothermic lithium 6 capture adds 4.78 MeV per bred triton.",
    ),
)

_ITER = MachineCase(
    name="ITER",
    description=(
        "The conventional low-field reference: a superconducting tokamak sized for a "
        "fusion gain of ten in inductively driven ELMy H-mode operation."
    ),
    source=(
        "ITER Physics Basis, Nuclear Fusion 39 (1999) 2137, and ITER Physics Basis "
        "Chapter 2, Nuclear Fusion 39 (1999) 2175. The inductive Q equal to ten "
        "baseline scenario."
    ),
    state=PlasmaState(
        geometry=PlasmaGeometry(
            major_radius=6.2, minor_radius=2.0, elongation=1.70, triangularity=0.33
        ),
        toroidal_field=5.3,
        plasma_current_ma=15.0,
        electron_density=1.0e20,
        temperature_kev=8.8,
        composition=PlasmaComposition(
            helium_fraction=0.04,
            impurities=(
                ImpurityRadiator(
                    species="beryllium",
                    atomic_number=4,
                    concentration=0.02,
                    loss_parameter=CORONAL_LOSS_PARAMETER["beryllium"],
                ),
                ImpurityRadiator(
                    species="tungsten",
                    atomic_number=74,
                    concentration=1.0e-5,
                    loss_parameter=CORONAL_LOSS_PARAMETER["tungsten"],
                ),
            ),
            mass_number=2.5,
        ),
        confinement_multiplier=1.0,
        wall_reflectivity=0.9,
    ),
    published=PublishedValues(
        fusion_power_mw=500.0,
        fusion_gain=10.0,
        confinement_time_s=3.7,
        normalised_beta=1.8,
        bootstrap_fraction=0.15,
        greenwald_fraction=0.85,
        auxiliary_power_mw=50.0,
    ),
    plant=None,
    notes=(
        "ITER has no power conversion system, so no electrical accounting is produced "
        "for it. Its value here is as the point the IPB98(y,2) scaling was constructed "
        "to project, which makes it the natural check on the confinement model.",
        "The published confinement time of 3.7 s corresponds to a loss power near "
        "87 MW. This model solves for the loss power rather than assuming it, so the "
        "two are compared without either being imposed.",
    ),
)

_SPARC = MachineCase(
    name="SPARC",
    description=(
        "A compact high-field deuterium tritium tokamak sized to demonstrate a fusion "
        "gain above one, and the experimental step preceding ARC."
    ),
    source=(
        "A. J. Creely et al., 'Overview of the SPARC tokamak', Journal of Plasma "
        "Physics 86 (2020) 865860502."
    ),
    state=PlasmaState(
        geometry=PlasmaGeometry(
            major_radius=1.85, minor_radius=0.57, elongation=1.97, triangularity=0.54
        ),
        toroidal_field=12.2,
        plasma_current_ma=8.7,
        electron_density=3.1e20,
        temperature_kev=7.3,
        composition=PlasmaComposition(helium_fraction=0.03, mass_number=2.5),
        confinement_multiplier=1.0,
        wall_reflectivity=0.9,
    ),
    published=PublishedValues(
        fusion_power_mw=140.0,
        fusion_gain=11.0,
        confinement_time_s=0.77,
        greenwald_fraction=0.37,
        auxiliary_power_mw=11.0,
    ),
    plant=None,
    notes=(
        "SPARC is a pulsed experiment with no breeding blanket and no power conversion, "
        "so only the plasma quantities are compared.",
        "It is the highest field point in this set, at 12.2 T, and therefore the point "
        "where the synchrotron term of this model carries the most weight.",
    ),
)

MACHINES: Final[dict[str, MachineCase]] = {
    case.name: case for case in (_ARC, _ITER, _SPARC)
}
"""Every published design point, keyed by machine name."""


def machine(name: str) -> MachineCase:
    """Return one design point by name.

    Raises:
        KeyError: If the name is not one of :data:`MACHINES`.
    """
    if name not in MACHINES:
        raise KeyError(f"unknown machine {name!r}; have {sorted(MACHINES)}")
    return MACHINES[name]
