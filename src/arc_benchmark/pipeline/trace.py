"""Structured records produced by the pipeline layer.

Everything the pipeline computes lands in one of these objects. They are frozen,
they hold no behaviour beyond derived properties, and they carry the inputs
alongside the outputs so that a report or a figure never has to reconstruct the
configuration that produced a number.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from arc_benchmark.algorithm.constraints import ConstraintReport
from arc_benchmark.algorithm.operating import OperatingPoint
from arc_benchmark.model.plant import PlantResult

__all__ = [
    "BenchmarkResult",
    "BenchmarkRow",
    "BenchmarkTrace",
    "SweepPoint",
    "SweepTrace",
]


@dataclass(frozen=True, slots=True)
class SweepPoint:
    """One evaluated point of a parametric sweep.

    Attributes:
        value: The swept independent variable at this point.
        point: The solved operating point.
        constraints: Every operational limit evaluated there.
        plant: The electrical accounting, or ``None`` when no plant model was
            supplied.
    """

    value: float
    point: OperatingPoint
    constraints: ConstraintReport
    plant: PlantResult | None

    @property
    def feasible(self) -> bool:
        """True when every operational limit is respected at this point."""
        return self.constraints.satisfied

    @property
    def binding_constraint(self) -> str:
        """Name of the constraint closest to its boundary here."""
        return self.constraints.binding.name


@dataclass(frozen=True, slots=True)
class SweepTrace:
    """A one-dimensional parametric sweep.

    Attributes:
        variable: Name of the swept variable.
        units: Units of the swept variable.
        policy: What was held fixed while it was swept. Stated explicitly
            because the answer to a sweep is meaningless without it.
        scaling_name: Confinement scaling used throughout.
        points: The evaluated points, in the order they were swept.
    """

    variable: str
    units: str
    policy: str
    scaling_name: str
    points: tuple[SweepPoint, ...]

    @property
    def values(self) -> tuple[float, ...]:
        """The swept independent variable at each point."""
        return tuple(p.value for p in self.points)

    @property
    def feasible_points(self) -> tuple[SweepPoint, ...]:
        """Only the points that respect every operational limit."""
        return tuple(p for p in self.points if p.feasible)

    def first_infeasible(self) -> SweepPoint | None:
        """The first point in sweep order that violates a limit, if any."""
        for point in self.points:
            if not point.feasible:
                return point
        return None


@dataclass(frozen=True, slots=True)
class BenchmarkRow:
    """One computed quantity compared against one published value.

    Attributes:
        quantity: What is being compared, in words.
        units: Units of both values.
        computed: What this model produces.
        published: The published value.
    """

    quantity: str
    units: str
    computed: float
    published: float

    @property
    def absolute_error(self) -> float:
        """Computed less published, in the units of the row."""
        return self.computed - self.published

    @property
    def relative_error(self) -> float:
        """Computed less published, divided by published.

        Infinite when the published value is zero, which no row here uses.
        """
        if self.published == 0.0:
            return math.inf
        return self.absolute_error / self.published


@dataclass(frozen=True, slots=True)
class BenchmarkResult:
    """One published design point reproduced, and the discrepancies found.

    Attributes:
        machine: Name of the design point.
        source: Citation for the published values.
        point: The solved operating point.
        constraints: Operational limits evaluated at that point.
        plant: Electrical accounting, or ``None`` for a device with no power
            conversion system.
        rows: One row per compared quantity.
        implied_confinement_multiplier: The H factor the confinement scaling
            would need in order to reproduce the published fusion power and
            auxiliary power together, given this model's stored energy and
            radiated power. ``None`` when the source does not quote both, or when
            the implied loss power is not positive. Comparing this against the H
            factor the source assumes separates a confinement disagreement from
            a fusion power disagreement, which the gain alone cannot do.
        notes: Statements about where and why this model is expected to differ,
            attached to the result rather than left to the reader.
    """

    machine: str
    source: str
    point: OperatingPoint
    constraints: ConstraintReport
    plant: PlantResult | None
    rows: tuple[BenchmarkRow, ...]
    implied_confinement_multiplier: float | None
    notes: tuple[str, ...]

    @property
    def worst_row(self) -> BenchmarkRow:
        """The row with the largest relative discrepancy."""
        return max(self.rows, key=lambda r: abs(r.relative_error))


@dataclass(frozen=True, slots=True)
class BenchmarkTrace:
    """Every benchmark case in one run.

    Attributes:
        scaling_name: Confinement scaling used.
        results: One result per design point.
    """

    scaling_name: str
    results: tuple[BenchmarkResult, ...]

    def named(self, machine: str) -> BenchmarkResult:
        """Return the result for one machine.

        Raises:
            KeyError: If no such machine was benchmarked.
        """
        for result in self.results:
            if result.machine == machine:
                return result
        raise KeyError(f"no benchmark for {machine!r}; have {[r.machine for r in self.results]}")
