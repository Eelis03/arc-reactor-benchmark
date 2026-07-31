"""Pipeline: named design points, parametric sweeps, and benchmark cases.

This layer wires the model and the algorithm together and produces a structured
trace. It does no plotting and prints nothing.
"""

from __future__ import annotations

from arc_benchmark.pipeline.benchmark import run_benchmark, run_benchmark_case
from arc_benchmark.pipeline.machines import MACHINES, MachineCase, PublishedValues, machine
from arc_benchmark.pipeline.sweep import (
    SweepInvariant,
    density_sweep,
    field_sweep,
    radius_sweep,
)
from arc_benchmark.pipeline.trace import (
    BenchmarkResult,
    BenchmarkRow,
    BenchmarkTrace,
    SweepPoint,
    SweepTrace,
)

__all__ = [
    "MACHINES",
    "BenchmarkResult",
    "BenchmarkRow",
    "BenchmarkTrace",
    "MachineCase",
    "PublishedValues",
    "SweepInvariant",
    "SweepPoint",
    "SweepTrace",
    "density_sweep",
    "field_sweep",
    "machine",
    "radius_sweep",
    "run_benchmark",
    "run_benchmark_case",
]
