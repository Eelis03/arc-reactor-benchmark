"""Figures. The only module in the package that touches matplotlib or the disk.

The non-interactive Agg backend is selected on import so that the examples run
identically with and without a display, which is what makes the integration tests
runnable in continuous integration.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import matplotlib
import numpy as np

from arc_benchmark.algorithm.balance import PlasmaComposition
from arc_benchmark.algorithm.lawson import lawson_triple_product
from arc_benchmark.algorithm.operating import OperatingPoint
from arc_benchmark.analysis.metrics import sweep_series
from arc_benchmark.pipeline.trace import SweepTrace

matplotlib.use("Agg")

from matplotlib import pyplot as plt

__all__ = ["plot_lawson_diagram", "plot_power_balance", "plot_sweep"]


def _save(figure: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(path, dpi=140)
    plt.close(figure)
    return path


def plot_sweep(trace: SweepTrace, path: Path) -> Path:
    """Plot fusion gain and net efficiency against the swept variable.

    Infeasible points are marked so that the reader cannot read a gain off a
    design point that violates a limit.

    Args:
        trace: The sweep to plot.
        path: Where to write the figure.

    Returns:
        The path written.
    """
    values = np.asarray(trace.values, dtype=np.float64)
    gains = np.asarray(sweep_series(trace, "gain"), dtype=np.float64)
    fusion = np.asarray(sweep_series(trace, "fusion_power"), dtype=np.float64)
    feasible = np.asarray([point.feasible for point in trace.points], dtype=bool)
    has_plant = all(point.plant is not None for point in trace.points)

    rows = 3 if has_plant else 2
    figure, axes = plt.subplots(rows, 1, figsize=(8.0, 3.0 * rows), sharex=True)

    axes[0].plot(values, fusion, color="#1f4e79")
    axes[0].set_ylabel("fusion power (MW)")
    axes[0].set_yscale("log")

    finite = np.isfinite(gains)
    axes[1].plot(values[finite], gains[finite], color="#8c2d04")
    axes[1].set_ylabel("fusion gain Q")
    axes[1].set_yscale("log")

    if has_plant:
        efficiency = np.asarray(sweep_series(trace, "net_efficiency"), dtype=np.float64)
        axes[2].plot(values, efficiency, color="#22633b")
        axes[2].axhline(0.0, color="0.5", linewidth=0.8)
        axes[2].set_ylabel("net electric / fusion")

    for axis in axes:
        axis.grid(True, alpha=0.3)
        if np.any(~feasible):
            axis.fill_between(
                values,
                *axis.get_ylim(),
                where=~feasible,
                color="#b03030",
                alpha=0.12,
                step="mid",
            )

    axes[-1].set_xlabel(f"{trace.variable} ({trace.units})")
    axes[0].set_title(f"{trace.variable} sweep, {trace.scaling_name}\n{trace.policy}", fontsize=9)
    return _save(figure, path)


def plot_lawson_diagram(
    composition: PlasmaComposition,
    points: tuple[tuple[str, OperatingPoint], ...],
    path: Path,
    gains: tuple[float, ...] = (1.0, 10.0, float("inf")),
    temperature_points: int = 240,
) -> Path:
    """Plot the triple product requirement against temperature, with design points.

    Args:
        composition: Composition the requirement curves are computed for.
        points: Named operating points to mark on the diagram.
        path: Where to write the figure.
        gains: Fusion gains to draw a requirement curve for. Infinity is
            ignition.
        temperature_points: Resolution of the requirement curves.

    Returns:
        The path written.
    """
    temperatures = np.linspace(3.0, 60.0, temperature_points, dtype=np.float64)
    figure, axis = plt.subplots(figsize=(8.0, 5.0))

    for gain in gains:
        values = np.array(
            [
                lawson_triple_product(float(t), composition, gain).triple_product
                for t in temperatures
            ],
            dtype=np.float64,
        )
        finite = np.isfinite(values)
        label = "ignition" if gain == float("inf") else f"Q = {gain:g}"
        axis.plot(temperatures[finite], values[finite], label=label)

    for name, point in points:
        axis.scatter(
            [point.state.temperature_kev],
            [point.triple_product],
            marker="o",
            zorder=5,
        )
        axis.annotate(
            name,
            (point.state.temperature_kev, point.triple_product),
            textcoords="offset points",
            xytext=(6, 4),
            fontsize=9,
        )

    axis.set_yscale("log")
    axis.set_xlabel("temperature (keV)")
    axis.set_ylabel("triple product n T tau (m^-3 keV s)")
    axis.set_title("Lawson requirement and the benchmarked design points", fontsize=10)
    axis.grid(True, alpha=0.3)
    axis.legend()
    return _save(figure, path)


def plot_power_balance(
    labelled_points: tuple[tuple[str, OperatingPoint], ...],
    path: Path,
) -> Path:
    """Stacked bars of sources and sinks for several operating points.

    Args:
        labelled_points: Named operating points to draw.
        path: Where to write the figure.

    Returns:
        The path written.
    """
    names = [name for name, _ in labelled_points]
    positions = np.arange(len(names), dtype=np.float64)
    width = 0.38

    alpha = np.array([p.terms.alpha_power_mw for _, p in labelled_points])
    auxiliary = np.array([max(p.auxiliary_power_mw, 0.0) for _, p in labelled_points])
    brem = np.array([p.terms.bremsstrahlung_mw for _, p in labelled_points])
    syn = np.array([p.terms.synchrotron_mw for _, p in labelled_points])
    line = np.array([p.terms.line_radiation_mw for _, p in labelled_points])
    transport = np.array([p.transport_power_mw for _, p in labelled_points])

    figure, axis = plt.subplots(figsize=(8.0, 5.0))
    axis.bar(positions - width / 2, alpha, width, label="alpha heating")
    axis.bar(positions - width / 2, auxiliary, width, bottom=alpha, label="auxiliary heating")
    axis.bar(positions + width / 2, brem, width, label="bremsstrahlung")
    axis.bar(positions + width / 2, syn, width, bottom=brem, label="synchrotron")
    axis.bar(positions + width / 2, line, width, bottom=brem + syn, label="line radiation")
    axis.bar(
        positions + width / 2,
        transport,
        width,
        bottom=brem + syn + line,
        label="transport",
    )

    axis.set_xticks(positions)
    axis.set_xticklabels(names)
    axis.set_ylabel("power (MW)")
    axis.set_title("Sources on the left of each pair, sinks on the right", fontsize=10)
    axis.legend(fontsize=8)
    axis.grid(True, axis="y", alpha=0.3)
    return _save(figure, path)
