"""Quantities extracted from a trace: fitted exponents and discrepancies.

The exponent fit exists because the argument for the high-field pathway is an
exponent. Stating that fusion power rises steeply with field is not a result;
measuring the exponent, saying what was held fixed while it was measured, and
checking it against the analytic expectation is.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass

import numpy as np

from arc_benchmark.pipeline.trace import SweepTrace

__all__ = [
    "PowerLawFit",
    "binding_constraint_counts",
    "fit_power_law",
    "sweep_series",
]


@dataclass(frozen=True, slots=True)
class PowerLawFit:
    """A least-squares fit of ``y = c x**k`` in log-log space.

    Attributes:
        exponent: The fitted ``k``.
        coefficient: The fitted ``c``.
        r_squared: Coefficient of determination of the fit in log space. Exactly
            one, to rounding, when the data really is a power law.
        points: Number of points the fit used.
    """

    exponent: float
    coefficient: float
    r_squared: float
    points: int


def fit_power_law(x: Sequence[float], y: Sequence[float]) -> PowerLawFit:
    """Fit ``y = c x**k`` by least squares on the logarithms.

    Both series must be strictly positive, since the fit is taken in log space.
    A plain polynomial fit of degree one on ``(log x, log y)`` is used rather
    than a nonlinear fit in the original space: the log-space fit weights
    relative error uniformly, which is what an exponent is a statement about.

    Args:
        x: Independent variable, strictly positive.
        y: Dependent variable, strictly positive.

    Returns:
        The fit, with its coefficient of determination.

    Raises:
        ValueError: If the series differ in length, have fewer than two points,
            or contain a non-positive value.
    """
    xs = np.asarray(x, dtype=np.float64)
    ys = np.asarray(y, dtype=np.float64)
    if xs.shape != ys.shape:
        raise ValueError(f"series must have the same shape, got {xs.shape} and {ys.shape}")
    if xs.size < 2:
        raise ValueError(f"need at least two points to fit an exponent, got {xs.size}")
    if np.any(xs <= 0.0) or np.any(ys <= 0.0):
        raise ValueError("both series must be strictly positive to fit in log space")

    log_x = np.log(xs)
    log_y = np.log(ys)
    slope, intercept = np.polyfit(log_x, log_y, 1)

    predicted = slope * log_x + intercept
    residual = float(np.sum((log_y - predicted) ** 2))
    total = float(np.sum((log_y - np.mean(log_y)) ** 2))
    r_squared = 1.0 if total == 0.0 else 1.0 - residual / total

    return PowerLawFit(
        exponent=float(slope),
        coefficient=float(np.exp(intercept)),
        r_squared=r_squared,
        points=int(xs.size),
    )


def sweep_series(trace: SweepTrace, quantity: str) -> tuple[float, ...]:
    """Extract one named quantity from every point of a sweep.

    Args:
        trace: The sweep.
        quantity: One of ``"fusion_power"``, ``"gain"``, ``"auxiliary_power"``,
            ``"confinement_time"``, ``"radiated_power"``, ``"loss_power"``,
            ``"triple_product"``, ``"net_electric"``, ``"net_efficiency"``,
            ``"recirculating_fraction"``, or ``"engineering_gain"``.

    Returns:
        The series, in sweep order.

    Raises:
        KeyError: If the quantity is not one of the supported names.
        ValueError: If a plant quantity is requested from a sweep run without a
            plant model.
    """
    plasma = {
        "fusion_power": lambda p: p.point.terms.fusion_power_mw,
        "gain": lambda p: p.point.fusion_gain,
        "auxiliary_power": lambda p: p.point.auxiliary_power_mw,
        "confinement_time": lambda p: p.point.confinement_time_s,
        "radiated_power": lambda p: p.point.terms.radiated_power_mw,
        "loss_power": lambda p: p.point.loss_power_mw,
        "triple_product": lambda p: p.point.triple_product,
    }
    plant = {
        "net_electric": lambda r: r.net_electric_mw,
        "net_efficiency": lambda r: r.net_efficiency,
        "recirculating_fraction": lambda r: r.recirculating_fraction,
        "engineering_gain": lambda r: r.engineering_gain,
    }

    if quantity in plasma:
        return tuple(plasma[quantity](p) for p in trace.points)
    if quantity in plant:
        if any(p.plant is None for p in trace.points):
            raise ValueError(
                f"quantity {quantity!r} needs a plant model, which this sweep was run without"
            )
        return tuple(plant[quantity](p.plant) for p in trace.points if p.plant is not None)
    raise KeyError(
        f"unknown quantity {quantity!r}; have {sorted(plasma) + sorted(plant)}"
    )


def binding_constraint_counts(trace: SweepTrace) -> dict[str, int]:
    """How many points of a sweep each constraint is the binding one at.

    Ordered by count, descending, then by name so the result is deterministic.
    """
    counts: dict[str, int] = {}
    for point in trace.points:
        name = point.binding_constraint
        counts[name] = counts.get(name, 0) + 1
    return dict(sorted(counts.items(), key=lambda item: (-item[1], item[0])))
