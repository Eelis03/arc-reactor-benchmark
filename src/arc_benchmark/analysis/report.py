"""Text tables built from traces.

Every function here returns a list of lines rather than printing. The examples
print them; the tests read them. Nothing in this module writes to a stream.
"""

from __future__ import annotations

import math

from arc_benchmark.algorithm.constraints import ConstraintReport
from arc_benchmark.algorithm.operating import OperatingPoint
from arc_benchmark.analysis.metrics import binding_constraint_counts
from arc_benchmark.pipeline.trace import BenchmarkResult, BenchmarkTrace, SweepTrace

__all__ = [
    "benchmark_lines",
    "constraint_lines",
    "power_balance_lines",
    "sweep_lines",
]


def _number(value: float, digits: int = 3) -> str:
    """Format a float for a table, keeping infinities readable."""
    if math.isinf(value):
        return "ignited" if value > 0.0 else "-inf"
    if math.isnan(value):
        return "n/a"
    return f"{value:.{digits}f}"


def power_balance_lines(point: OperatingPoint) -> list[str]:
    """A readable audit of one solved operating point."""
    terms = point.terms
    lines = [
        f"Operating point solved with {point.scaling_name}, "
        f"{point.convention.value} loss power convention",
        f"  volume                       {point.state.geometry.volume:10.2f} m^3",
        f"  electron density             {point.state.electron_density:10.3e} m^-3",
        f"  temperature                  {point.state.temperature_kev:10.2f} keV",
        f"  reactivity                   {terms.reactivity_m3_s:10.4e} m^3/s",
        "  sources",
        f"    alpha heating              {terms.alpha_power_mw:10.2f} MW",
        f"    auxiliary heating          {point.auxiliary_power_mw:10.2f} MW",
        "  sinks",
        f"    bremsstrahlung             {terms.bremsstrahlung_mw:10.2f} MW",
        f"    synchrotron                {terms.synchrotron_mw:10.2f} MW",
        f"    line radiation             {terms.line_radiation_mw:10.2f} MW",
        f"    transport                  {point.transport_power_mw:10.2f} MW",
        f"  closure residual             {point.residual_mw:10.3e} MW",
        f"  fusion power                 {terms.fusion_power_mw:10.2f} MW",
        f"  neutron power                {terms.neutron_power_mw:10.2f} MW",
        f"  stored energy                {terms.stored_energy_mj:10.2f} MJ",
        f"  confinement time             {point.confinement_time_s:10.3f} s",
        f"  loss power                   {point.loss_power_mw:10.2f} MW",
        f"  fusion gain Q                {_number(point.fusion_gain, 2):>10}",
        f"  radiated fraction            {point.radiated_fraction:10.3f}",
        f"  triple product               {point.triple_product:10.4e} m^-3 keV s",
        f"  ignited                      {point.ignited!s:>10}",
    ]
    return lines


def constraint_lines(report: ConstraintReport) -> list[str]:
    """One line per constraint, with the verdict and the margin."""
    lines = [
        f"{'constraint':16} {'value':>12} {'limit':>12} {'used':>8} verdict",
    ]
    for check in report.checks:
        verdict = "ok" if check.satisfied else "VIOLATED"
        lines.append(
            f"{check.name:16} {check.value:12.4g} {check.limit:12.4g} "
            f"{check.utilisation:8.3f} {verdict}"
        )
    lines.append(
        f"binding constraint: {report.binding.name} at "
        f"{report.binding.utilisation:.3f} of its limit"
    )
    if report.violations:
        names = ", ".join(check.name for check in report.violations)
        lines.append(f"VIOLATED: {names}")
    else:
        lines.append("all constraints satisfied")
    return lines


def _benchmark_result_lines(result: BenchmarkResult) -> list[str]:
    lines = [
        f"{result.machine}",
        f"  source: {result.source}",
        f"  {'quantity':26} {'units':>12} {'computed':>12} {'published':>12} {'error':>9}",
    ]
    for row in result.rows:
        lines.append(
            f"  {row.quantity:26} {row.units:>12} {row.computed:12.3f} "
            f"{row.published:12.3f} {row.relative_error * 100:+8.1f}%"
        )
    assumed = result.point.state.confinement_multiplier
    implied = result.implied_confinement_multiplier
    if implied is not None:
        lines.append(
            f"  H factor assumed by the source {assumed:.3f}, implied by the published "
            f"operating point {implied:.3f}, ratio {implied / assumed:.3f}"
        )
    lines.append(
        f"  constraints: {'all satisfied' if result.constraints.satisfied else 'VIOLATED'}"
        f", binding {result.constraints.binding.name} at "
        f"{result.constraints.binding.utilisation:.3f}"
    )
    for note in result.notes:
        lines.append(f"  note: {note}")
    return lines


def benchmark_lines(trace: BenchmarkTrace) -> list[str]:
    """A comparison table for every benchmarked design point."""
    lines = [f"Benchmark against published design points, scaling {trace.scaling_name}"]
    for result in trace.results:
        lines.append("")
        lines.extend(_benchmark_result_lines(result))
    return lines


def sweep_lines(trace: SweepTrace, stride: int = 1) -> list[str]:
    """A table of one sweep, one row per point.

    Args:
        trace: The sweep.
        stride: Print every ``stride`` th point. A fine sweep is useful for
            fitting an exponent and unreadable as a table, so the table is
            thinned rather than the sweep.

    Returns:
        The table as lines.

    Raises:
        ValueError: If the stride is not positive.
    """
    if stride < 1:
        raise ValueError(f"stride must be at least 1, got {stride}")

    has_plant = all(point.plant is not None for point in trace.points)
    header = f"{trace.variable} sweep, {trace.units}, {trace.scaling_name}, holding: {trace.policy}"
    columns = (
        f"{'value':>10} {'P_fus MW':>10} {'P_aux MW':>10} {'Q':>10} {'tau_E s':>9} "
        f"{'binding':>14} {'ok':>5}"
    )
    if has_plant:
        columns += f" {'net MWe':>9} {'eta_net':>8}"
    lines = [header, columns]

    for point in trace.points[::stride]:
        solved = point.point
        row = (
            f"{point.value:10.4g} {solved.terms.fusion_power_mw:10.2f} "
            f"{solved.auxiliary_power_mw:10.2f} {_number(solved.fusion_gain, 2):>10} "
            f"{solved.confinement_time_s:9.3f} {point.binding_constraint:>14} "
            f"{('yes' if point.feasible else 'NO'):>5}"
        )
        if has_plant and point.plant is not None:
            row += f" {point.plant.net_electric_mw:9.2f} {point.plant.net_efficiency:8.4f}"
        lines.append(row)

    counts = binding_constraint_counts(trace)
    summary = ", ".join(f"{name} at {count} points" for name, count in counts.items())
    lines.append(f"binding constraint across the sweep: {summary}")

    # The lowest and the highest feasible value, said in those words rather than
    # as a band, because nothing here guarantees that the feasible points are
    # contiguous and a report should not imply that they are.
    feasible = trace.feasible_points
    if feasible:
        lines.append(
            f"feasible at {len(feasible)} of {len(trace.points)} points, lowest "
            f"{feasible[0].value:.4g} {trace.units}, highest {feasible[-1].value:.4g} "
            f"{trace.units}"
        )
    else:
        lines.append("no point in the sweep satisfies every constraint")

    first_bad = trace.first_infeasible()
    if first_bad is None:
        lines.append("every point in the sweep satisfies every constraint")
    else:
        names = ", ".join(check.name for check in first_bad.constraints.violations)
        lines.append(
            f"first infeasible point at {trace.variable} = {first_bad.value:.4g} "
            f"{trace.units}, violating {names}"
        )
    return lines
