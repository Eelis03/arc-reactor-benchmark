"""Reproduction of published design points, with the discrepancies reported.

The point of a benchmark is to find where a model is wrong, so every quantity the
source quotes and this model computes is compared, and the comparison is
reported whether it agrees or not. A benchmark that only reported agreements
would be a decoration.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass

from arc_benchmark.algorithm.balance import LossPowerConvention
from arc_benchmark.algorithm.constraints import ConstraintLimits, evaluate_constraints
from arc_benchmark.algorithm.operating import OperatingPoint, solve_operating_point
from arc_benchmark.algorithm.protocols import ConfinementScaling
from arc_benchmark.model.confinement import IPB98Y2
from arc_benchmark.model.constants import ALPHA_FRACTION
from arc_benchmark.model.limits import greenwald_density
from arc_benchmark.model.plant import PlantResult, evaluate_plant
from arc_benchmark.model.profiles import ProfileShape
from arc_benchmark.pipeline.machines import MACHINES, MachineCase, PublishedValues
from arc_benchmark.pipeline.trace import BenchmarkResult, BenchmarkRow, BenchmarkTrace

__all__ = [
    "evaluate_plant_for_point",
    "implied_confinement_multiplier",
    "run_benchmark",
    "run_benchmark_case",
]


@dataclass(frozen=True, slots=True)
class _Comparison:
    """Internal description of one comparable quantity."""

    quantity: str
    units: str
    computed: float
    published: float | None


def evaluate_plant_for_point(
    point: OperatingPoint,
    case: MachineCase,
) -> PlantResult | None:
    """Run the electrical accounting for a solved point, if the case has a plant.

    The auxiliary power passed to the plant model is floored at zero. A point
    that solves to a negative auxiliary power is past ignition, and negative
    heating draw is not a thing a plant can do. The floor is a modelling
    limitation rather than a physical statement: a steady-state device still
    needs current drive power whether or not it needs heating power, and this
    model does not separate the two.
    """
    if case.plant is None:
        return None
    return evaluate_plant(
        fusion_power_mw=point.terms.fusion_power_mw,
        alpha_power_mw=point.terms.alpha_power_mw,
        neutron_power_mw=point.terms.neutron_power_mw,
        auxiliary_power_mw=max(point.auxiliary_power_mw, 0.0),
        parameters=case.plant,
    )


def implied_confinement_multiplier(
    point: OperatingPoint,
    published: PublishedValues,
    scaling: ConfinementScaling,
) -> float | None:
    """H factor needed to reproduce a published operating point.

    The published fusion power fixes the alpha heating, the published auxiliary
    power fixes the external heating, and this model supplies the radiated power
    and the stored energy. Those together fix the loss power and therefore the
    confinement time the design point requires. Dividing that by what the scaling
    predicts at the same loss power gives the H factor the point implies.

    This separates two failures that the fusion gain alone conflates. If the
    implied H factor matches the one the source assumes, then the model's
    confinement physics agrees with the source and any disagreement in gain comes
    from the fusion power. If it does not, the disagreement is in the confinement
    or in the radiated power.

    Args:
        point: The solved operating point, used for its stored energy, its
            radiated power, and its geometry.
        published: The design values.
        scaling: The confinement scaling to measure against.

    Returns:
        The implied H factor, or ``None`` when the source does not quote both a
        fusion power and an auxiliary power, or when the implied loss power is
        not positive.
    """
    if published.fusion_power_mw is None or published.auxiliary_power_mw is None:
        return None

    alpha_mw = ALPHA_FRACTION * published.fusion_power_mw
    radiated = (
        point.terms.radiated_power_mw if point.convention is LossPowerConvention.SEPARATRIX else 0.0
    )
    loss_power_mw = alpha_mw + published.auxiliary_power_mw - radiated
    if loss_power_mw <= 0.0:
        return None

    required_tau = point.terms.stored_energy_mj / loss_power_mw
    predicted_tau = scaling.tau_e(point.state.confinement_inputs(loss_power_mw))
    return required_tau / predicted_tau


def _rows(
    point: OperatingPoint,
    plant: PlantResult | None,
    published: PublishedValues,
    greenwald_fraction: float,
    beta_n: float,
    bootstrap: float,
) -> tuple[BenchmarkRow, ...]:
    comparisons = (
        _Comparison("fusion power", "MW", point.terms.fusion_power_mw, published.fusion_power_mw),
        _Comparison("fusion gain Q", "", point.fusion_gain, published.fusion_gain),
        _Comparison(
            "energy confinement time",
            "s",
            point.confinement_time_s,
            published.confinement_time_s,
        ),
        _Comparison(
            "auxiliary power",
            "MW",
            point.auxiliary_power_mw,
            published.auxiliary_power_mw,
        ),
        _Comparison("normalised beta", "% m T / MA", beta_n, published.normalised_beta),
        _Comparison("bootstrap fraction", "", bootstrap, published.bootstrap_fraction),
        _Comparison("Greenwald fraction", "", greenwald_fraction, published.greenwald_fraction),
        _Comparison(
            "thermal power",
            "MW",
            plant.thermal_power_mw if plant is not None else math.nan,
            published.thermal_power_mw,
        ),
        _Comparison(
            "net electric power",
            "MWe",
            plant.net_electric_mw if plant is not None else math.nan,
            published.net_electric_mw,
        ),
    )
    return tuple(
        BenchmarkRow(
            quantity=c.quantity,
            units=c.units,
            computed=c.computed,
            published=c.published,
        )
        for c in comparisons
        if c.published is not None and math.isfinite(c.computed)
    )


def run_benchmark_case(
    case: MachineCase,
    scaling: ConfinementScaling = IPB98Y2,
    convention: LossPowerConvention = LossPowerConvention.SEPARATRIX,
    limits: ConstraintLimits | None = None,
    profile: ProfileShape | None = None,
) -> BenchmarkResult:
    """Solve one published design point and compare it against its source.

    Args:
        case: The design point.
        scaling: Confinement scaling to solve with.
        convention: Loss power convention.
        limits: Constraint thresholds. Defaults to the published limits.
        profile: Profile shape to apply, overriding whatever the case carries.
            ``None`` leaves the case alone, which for every case defined here
            means flat profiles and a genuinely zero-dimensional answer.

    Returns:
        The solved point, its constraint report, its electrical accounting, and
        one comparison row per quantity the source quotes.
    """
    state = case.state if profile is None else dataclasses.replace(case.state, profile=profile)
    point = solve_operating_point(state, scaling, convention)
    report = evaluate_constraints(point, limits)
    plant = evaluate_plant_for_point(point, case)

    density_limit = greenwald_density(state.plasma_current_ma, state.geometry.minor_radius)
    # The sources quote the Greenwald fraction against the line-averaged density,
    # so the comparison is made against the same average.
    greenwald_fraction = state.line_averaged_density / density_limit
    beta_n = report.named("troyon").value
    bootstrap = report.named("bootstrap").value

    return BenchmarkResult(
        machine=case.name,
        source=case.source,
        point=point,
        constraints=report,
        plant=plant,
        rows=_rows(point, plant, case.published, greenwald_fraction, beta_n, bootstrap),
        implied_confinement_multiplier=implied_confinement_multiplier(
            point, case.published, scaling
        ),
        notes=case.notes,
    )


def run_benchmark(
    scaling: ConfinementScaling = IPB98Y2,
    convention: LossPowerConvention = LossPowerConvention.SEPARATRIX,
    machines: tuple[str, ...] | None = None,
    limits: ConstraintLimits | None = None,
    profile: ProfileShape | None = None,
) -> BenchmarkTrace:
    """Solve every published design point and compare each against its source.

    Args:
        scaling: Confinement scaling to solve with.
        convention: Loss power convention.
        machines: Names of the design points to run. Defaults to all of them, in
            the order they are declared.
        limits: Constraint thresholds.
        profile: Profile shape applied to every case, or ``None`` to leave each
            case as declared.

    Returns:
        One result per machine.

    Raises:
        KeyError: If a requested machine is not defined.
    """
    names = machines if machines is not None else tuple(MACHINES)
    results = []
    for name in names:
        if name not in MACHINES:
            raise KeyError(f"unknown machine {name!r}; have {sorted(MACHINES)}")
        results.append(run_benchmark_case(MACHINES[name], scaling, convention, limits, profile))
    return BenchmarkTrace(scaling_name=scaling.name, results=tuple(results))
