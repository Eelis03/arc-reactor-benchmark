"""Analysis: fitted exponents, discrepancy tables, and figures.

This layer reads what the pipeline produces. It never solves anything and never
reaches back into the model or the algorithm layers for physics.
"""

from __future__ import annotations

from arc_benchmark.analysis.metrics import (
    PowerLawFit,
    binding_constraint_counts,
    fit_power_law,
    sweep_series,
)
from arc_benchmark.analysis.report import (
    benchmark_lines,
    constraint_lines,
    power_balance_lines,
    sweep_lines,
)
from arc_benchmark.analysis.tables import benchmark_frame, constraint_frame, sweep_frame

__all__ = [
    "PowerLawFit",
    "benchmark_frame",
    "benchmark_lines",
    "binding_constraint_counts",
    "constraint_frame",
    "constraint_lines",
    "fit_power_law",
    "power_balance_lines",
    "sweep_frame",
    "sweep_lines",
    "sweep_series",
]
