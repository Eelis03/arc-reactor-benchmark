"""Tier one: constraint verdicts, sweeps, and the benchmark machinery.

The property this file exists to establish is that a design point violating a
published limit is reported as violating it, rather than returned as though it
were feasible.
"""

from __future__ import annotations

import dataclasses

import numpy as np
import pytest

from arc_benchmark.algorithm.constraints import (
    ConstraintLimits,
    ConstraintSense,
    cylindrical_q,
    evaluate_constraints,
)
from arc_benchmark.algorithm.operating import solve_operating_point
from arc_benchmark.analysis.metrics import fit_power_law, sweep_series
from arc_benchmark.analysis.report import (
    benchmark_lines,
    constraint_lines,
    power_balance_lines,
    sweep_lines,
)
from arc_benchmark.analysis.tables import benchmark_frame, constraint_frame, sweep_frame
from arc_benchmark.model.confinement import IPB98Y2, ITER89P
from arc_benchmark.model.limits import greenwald_density
from arc_benchmark.pipeline.benchmark import run_benchmark, run_benchmark_case
from arc_benchmark.pipeline.machines import MACHINES, machine
from arc_benchmark.pipeline.sweep import (
    SweepInvariant,
    density_sweep,
    field_sweep,
    radius_sweep,
)

_ARC = machine("ARC")
_ITER = machine("ITER")


def test_a_point_above_the_greenwald_limit_is_flagged() -> None:
    """Pushing the density past the limit produces a violation, not a silent answer."""
    limit = greenwald_density(_ITER.state.plasma_current_ma, _ITER.state.geometry.minor_radius)
    over = dataclasses.replace(_ITER.state, electron_density=1.2 * limit)
    report = evaluate_constraints(solve_operating_point(over, IPB98Y2))

    assert not report.satisfied
    greenwald = report.named("greenwald")
    assert not greenwald.satisfied
    assert greenwald.utilisation == pytest.approx(1.2, rel=1.0e-12)
    assert greenwald.margin < 0.0
    assert "greenwald" in {check.name for check in report.violations}


def test_a_point_below_the_safety_factor_limit_is_flagged() -> None:
    """A current large enough to drive the safety factor under two is a violation."""
    over_current = dataclasses.replace(_ITER.state, plasma_current_ma=30.0)
    report = evaluate_constraints(solve_operating_point(over_current, IPB98Y2))
    check = report.named("safety_factor")
    assert not check.satisfied
    assert check.sense is ConstraintSense.LOWER
    assert check.value < check.limit


def test_utilisation_reads_the_same_way_for_both_senses() -> None:
    """Above one is a violation whether the limit is an upper or a lower bound."""
    report = evaluate_constraints(solve_operating_point(_ARC.state, IPB98Y2))
    for check in report.checks:
        assert check.satisfied == (check.utilisation <= 1.0)
        assert check.margin == pytest.approx(1.0 - check.utilisation, rel=1.0e-15)


def test_the_binding_constraint_is_the_one_closest_to_its_boundary() -> None:
    """The report names the constraint with the highest utilisation."""
    report = evaluate_constraints(solve_operating_point(_ARC.state, IPB98Y2))
    highest = max(check.utilisation for check in report.checks)
    assert report.binding.utilisation == highest


def test_the_lh_threshold_check_is_dropped_for_an_l_mode_scaling() -> None:
    """An L-mode solve is not judged against an H-mode access threshold."""
    point = solve_operating_point(_ARC.state, ITER89P)
    with_check = evaluate_constraints(point, ConstraintLimits(require_h_mode=True))
    without_check = evaluate_constraints(point, ConstraintLimits(require_h_mode=False))
    assert "lh_threshold" in {check.name for check in with_check.checks}
    assert "lh_threshold" not in {check.name for check in without_check.checks}


def test_asking_for_an_unknown_constraint_raises() -> None:
    """Looking up a check that was not evaluated is an error, not a silent default."""
    report = evaluate_constraints(solve_operating_point(_ARC.state, IPB98Y2))
    with pytest.raises(KeyError, match="no constraint named"):
        report.named("not_a_constraint")


def test_cylindrical_safety_factor_is_reported_alongside_q95() -> None:
    """Both safety factors are available and the edge value is the larger."""
    point = solve_operating_point(_ARC.state, IPB98Y2)
    report = evaluate_constraints(point)
    assert cylindrical_q(point) > 0.0
    assert report.named("safety_factor").value > cylindrical_q(point)


@pytest.mark.parametrize(
    ("invariant", "expected_exponent"),
    [
        (SweepInvariant.FIXED_DENSITY, 0.0),
        (SweepInvariant.FIXED_GREENWALD_FRACTION, 2.0),
        (SweepInvariant.FIXED_BETA, 4.0),
    ],
)
def test_field_sweep_reproduces_the_analytic_exponent(
    invariant: SweepInvariant, expected_exponent: float
) -> None:
    """Each invariant produces exactly the field exponent it is derived to produce.

    Tolerance: the swept quantity is an exact power law in the field, so a
    least-squares fit on the logarithms recovers the exponent up to the
    conditioning of the normal equations, which for eleven well spread points is
    of order 1e-13. The tolerance is set at 1e-9, four orders above that, and is
    derived from the conditioning rather than from the residual observed.
    """
    fields = np.linspace(6.0, 12.0, 11)
    trace = field_sweep(_ARC.state, fields, invariant, IPB98Y2)
    fit = fit_power_law(trace.values, sweep_series(trace, "fusion_power"))
    assert fit.exponent == pytest.approx(expected_exponent, abs=1.0e-9)
    assert fit.r_squared == pytest.approx(1.0, abs=1.0e-12)


def test_radius_sweep_reproduces_its_analytic_exponent() -> None:
    """Fusion power goes as the first power of major radius at a fixed Greenwald fraction.

    The volume goes as the cube of the radius while the Greenwald density falls
    as its inverse, and fusion power carries the square of the density, so the
    exponent is ``3 - 2 = 1``.
    """
    radii = np.linspace(2.5, 6.5, 11)
    trace = radius_sweep(_ARC.state, radii, SweepInvariant.FIXED_GREENWALD_FRACTION, IPB98Y2)
    fit = fit_power_law(trace.values, sweep_series(trace, "fusion_power"))
    assert fit.exponent == pytest.approx(1.0, abs=1.0e-9)


def test_every_sweep_holds_the_safety_factor_it_says_it_holds() -> None:
    """The cylindrical safety factor is constant across a field or radius sweep.

    Tolerance: the current is solved from the safety factor and then the safety
    factor is recomputed from that current, a round trip of about ten
    floating-point operations, so 1e-12 relative is four orders above the
    rounding.
    """
    base_q = cylindrical_q(solve_operating_point(_ARC.state, IPB98Y2))
    for trace in (
        field_sweep(_ARC.state, np.linspace(6.0, 12.0, 9), scaling=IPB98Y2),
        radius_sweep(_ARC.state, np.linspace(2.5, 6.5, 9), scaling=IPB98Y2),
    ):
        for point in trace.points:
            assert cylindrical_q(point.point) == pytest.approx(base_q, rel=1.0e-12)


def test_density_sweep_changes_only_the_density() -> None:
    """Nothing else in the state moves while the density is swept."""
    densities = np.linspace(0.8e20, 2.0e20, 7)
    trace = density_sweep(_ARC.state, densities, IPB98Y2)
    for value, point in zip(densities, trace.points, strict=True):
        state = point.point.state
        assert state.electron_density == value
        assert state.toroidal_field == _ARC.state.toroidal_field
        assert state.plasma_current_ma == _ARC.state.plasma_current_ma
        assert state.temperature_kev == _ARC.state.temperature_kev


def test_a_sweep_reports_its_first_infeasible_point() -> None:
    """Pushing the density above the Greenwald limit makes the sweep say so."""
    limit = greenwald_density(_ARC.state.plasma_current_ma, _ARC.state.geometry.minor_radius)
    densities = np.linspace(0.5 * limit, 1.5 * limit, 21)
    trace = density_sweep(_ARC.state, densities, IPB98Y2)
    first_bad = trace.first_infeasible()
    assert first_bad is not None
    assert len(trace.feasible_points) < len(trace.points)


def test_sweep_carries_a_plant_result_only_when_a_plant_was_supplied() -> None:
    """A sweep without plant parameters produces no electrical accounting."""
    fields = np.linspace(8.0, 10.0, 3)
    without = field_sweep(_ARC.state, fields, scaling=IPB98Y2)
    with_plant = field_sweep(_ARC.state, fields, scaling=IPB98Y2, plant=_ARC.plant)
    assert all(point.plant is None for point in without.points)
    assert all(point.plant is not None for point in with_plant.points)
    with pytest.raises(ValueError, match="needs a plant model"):
        sweep_series(without, "net_efficiency")


def test_benchmark_produces_one_row_per_published_quantity() -> None:
    """Only quantities the source actually quotes are compared."""
    trace = run_benchmark(IPB98Y2)
    assert {result.machine for result in trace.results} == set(MACHINES)

    iter_result = trace.named("ITER")
    quantities = {row.quantity for row in iter_result.rows}
    assert "energy confinement time" in quantities
    assert "net electric power" not in quantities

    arc_result = trace.named("ARC")
    assert "net electric power" in {row.quantity for row in arc_result.rows}
    assert arc_result.plant is not None
    assert iter_result.plant is None


def test_benchmark_rows_report_signed_relative_error() -> None:
    """The reported error is computed minus published over published."""
    result = run_benchmark_case(_ARC, IPB98Y2)
    for row in result.rows:
        assert row.absolute_error == pytest.approx(row.computed - row.published, rel=1.0e-14)
        assert row.relative_error == pytest.approx(
            row.absolute_error / row.published, rel=1.0e-14
        )


def test_implied_confinement_multiplier_is_reported_for_every_machine() -> None:
    """Each source quotes both a fusion power and an auxiliary power, so all three resolve."""
    trace = run_benchmark(IPB98Y2)
    for result in trace.results:
        implied = result.implied_confinement_multiplier
        assert implied is not None
        assert 0.5 < implied < 5.0


def test_a_peaked_profile_raises_the_computed_fusion_power() -> None:
    """The profile override reaches the solve and changes the answer in the right direction."""
    from arc_benchmark.model.profiles import ProfileShape

    flat = run_benchmark_case(_ARC, IPB98Y2)
    peaked = run_benchmark_case(_ARC, IPB98Y2, profile=ProfileShape(0.4, 1.0))
    assert peaked.point.terms.fusion_power_mw > flat.point.terms.fusion_power_mw


def test_asking_for_an_unknown_machine_raises() -> None:
    """A misspelled machine name is an error, not an empty benchmark."""
    with pytest.raises(KeyError, match="unknown machine"):
        machine("TOKAMAK")
    with pytest.raises(KeyError, match="unknown machine"):
        run_benchmark(IPB98Y2, machines=("TOKAMAK",))


def test_power_law_fit_rejects_input_it_cannot_take_the_logarithm_of() -> None:
    """A non-positive value, or a single point, raises rather than returning a number."""
    with pytest.raises(ValueError, match="strictly positive"):
        fit_power_law([1.0, 2.0], [1.0, -1.0])
    with pytest.raises(ValueError, match="at least two points"):
        fit_power_law([1.0], [1.0])
    with pytest.raises(ValueError, match="same shape"):
        fit_power_law([1.0, 2.0], [1.0])


def test_reports_render_without_raising_and_mention_the_verdict() -> None:
    """Every report function produces non-empty text containing its key facts."""
    point = solve_operating_point(_ARC.state, IPB98Y2)
    report = evaluate_constraints(point)

    balance = power_balance_lines(point)
    assert any("closure residual" in line for line in balance)
    assert any("fusion gain" in line for line in balance)

    constraints = constraint_lines(report)
    assert any("binding constraint" in line for line in constraints)

    trace = field_sweep(
        _ARC.state, np.linspace(8.0, 12.0, 5), scaling=IPB98Y2, plant=_ARC.plant
    )
    sweep_text = sweep_lines(trace)
    assert any("binding constraint across the sweep" in line for line in sweep_text)
    assert any("eta_net" in line for line in sweep_text)

    benchmark_text = benchmark_lines(run_benchmark(IPB98Y2))
    assert any("ARC" in line for line in benchmark_text)
    assert any("H factor assumed by the source" in line for line in benchmark_text)


def test_frames_carry_the_same_numbers_as_the_trace() -> None:
    """The data frames are a view of the trace, not a rounded copy of it."""
    trace = field_sweep(
        _ARC.state, np.linspace(8.0, 12.0, 5), scaling=IPB98Y2, plant=_ARC.plant
    )
    frame = sweep_frame(trace)
    assert len(frame) == len(trace.points)
    assert list(frame["toroidal field on axis"]) == list(trace.values)
    for row, point in zip(frame.itertuples(), trace.points, strict=True):
        assert row.fusion_power_mw == point.point.terms.fusion_power_mw
        assert row.feasible == point.feasible
        assert row.binding_constraint == point.binding_constraint
        assert point.plant is not None
        assert row.net_efficiency == point.plant.net_efficiency

    constraints = constraint_frame(trace)
    assert len(constraints) == sum(len(p.constraints.checks) for p in trace.points)
    expected_names = {check.name for check in trace.points[0].constraints.checks}
    assert set(constraints["constraint"]) == expected_names

    benchmark = benchmark_frame(run_benchmark(IPB98Y2))
    assert set(benchmark["machine"]) == set(MACHINES)
    for row in benchmark.itertuples():
        assert row.absolute_error == pytest.approx(row.computed - row.published, rel=1.0e-14)


def test_a_sweep_without_a_plant_produces_a_frame_without_plant_columns() -> None:
    """The plant columns appear only when there is a plant result to put in them."""
    trace = field_sweep(_ARC.state, np.linspace(8.0, 12.0, 3), scaling=IPB98Y2)
    frame = sweep_frame(trace)
    assert "net_efficiency" not in frame.columns
    assert "fusion_power_mw" in frame.columns


def test_sweep_report_rejects_a_zero_stride() -> None:
    """A stride below one would produce an empty or reversed table and raises."""
    trace = field_sweep(_ARC.state, np.linspace(8.0, 12.0, 5), scaling=IPB98Y2)
    with pytest.raises(ValueError, match="stride"):
        sweep_lines(trace, stride=0)
