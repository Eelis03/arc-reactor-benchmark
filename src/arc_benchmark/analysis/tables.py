"""Traces rendered as data frames.

The text tables in :mod:`arc_benchmark.analysis.report` are for reading. These
are for working with: sorting a sweep by net efficiency, filtering it to the
feasible points, joining a benchmark against another study, or writing either to
a file. Nothing here reformats a number, so a value in a frame is the value the
solver produced and not a rounded copy of it.
"""

from __future__ import annotations

import pandas as pd

from arc_benchmark.pipeline.trace import BenchmarkTrace, SweepTrace

__all__ = ["benchmark_frame", "constraint_frame", "sweep_frame"]


def benchmark_frame(trace: BenchmarkTrace) -> pd.DataFrame:
    """One row per compared quantity, across every benchmarked machine.

    Columns: ``machine``, ``quantity``, ``units``, ``computed``, ``published``,
    ``absolute_error``, ``relative_error``.
    """
    records = [
        {
            "machine": result.machine,
            "quantity": row.quantity,
            "units": row.units,
            "computed": row.computed,
            "published": row.published,
            "absolute_error": row.absolute_error,
            "relative_error": row.relative_error,
        }
        for result in trace.results
        for row in result.rows
    ]
    return pd.DataFrame.from_records(
        records,
        columns=[
            "machine",
            "quantity",
            "units",
            "computed",
            "published",
            "absolute_error",
            "relative_error",
        ],
    )


def sweep_frame(trace: SweepTrace) -> pd.DataFrame:
    """One row per swept point, with the plasma solution and the plant result.

    The plant columns are present only when the sweep was run with plant
    parameters. The feasibility verdict and the binding constraint are carried
    alongside the numbers so that a filter on feasibility does not need a second
    pass over the trace.
    """
    records = []
    for point in trace.points:
        solved = point.point
        record: dict[str, object] = {
            trace.variable: point.value,
            "fusion_power_mw": solved.terms.fusion_power_mw,
            "alpha_power_mw": solved.terms.alpha_power_mw,
            "auxiliary_power_mw": solved.auxiliary_power_mw,
            "fusion_gain": solved.fusion_gain,
            "confinement_time_s": solved.confinement_time_s,
            "loss_power_mw": solved.loss_power_mw,
            "radiated_power_mw": solved.terms.radiated_power_mw,
            "triple_product": solved.triple_product,
            "feasible": point.feasible,
            "binding_constraint": point.binding_constraint,
        }
        if point.plant is not None:
            record.update(
                {
                    "thermal_power_mw": point.plant.thermal_power_mw,
                    "gross_electric_mw": point.plant.gross_electric_mw,
                    "recirculating_mw": point.plant.recirculating_mw,
                    "net_electric_mw": point.plant.net_electric_mw,
                    "net_efficiency": point.plant.net_efficiency,
                    "engineering_gain": point.plant.engineering_gain,
                }
            )
        records.append(record)
    return pd.DataFrame.from_records(records)


def constraint_frame(trace: SweepTrace) -> pd.DataFrame:
    """One row per swept point per constraint, in long form.

    Long form rather than wide because the set of constraints depends on whether
    the H-mode access threshold applies, and a wide frame would then change its
    columns between runs.
    """
    records = [
        {
            trace.variable: point.value,
            "constraint": check.name,
            "value": check.value,
            "limit": check.limit,
            "utilisation": check.utilisation,
            "satisfied": check.satisfied,
        }
        for point in trace.points
        for check in point.constraints.checks
    ]
    return pd.DataFrame.from_records(
        records,
        columns=[
            trace.variable,
            "constraint",
            "value",
            "limit",
            "utilisation",
            "satisfied",
        ],
    )
