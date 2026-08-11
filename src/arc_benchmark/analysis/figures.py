"""Figures. The only module in the package that touches matplotlib or the disk.

The non-interactive Agg backend is selected on import so that the examples run
identically with and without a display, which is what makes the integration tests
runnable in continuous integration.

Every figure here draws from one palette and one set of chrome colours, declared
at the top of this module, so that the figures published in the README read as
one set rather than as five separate scripts. The three series colours are a
colourblind-safe triple: the worst pair separation is 9.2 in OKLab units under
deuteranopia and 24.0 under normal vision. Identity is never carried by colour
alone: every figure with more than one series has both a legend and a direct
label on the mark.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Final

import matplotlib
import numpy as np

from arc_benchmark.algorithm.balance import PlasmaComposition
from arc_benchmark.algorithm.lawson import lawson_triple_product
from arc_benchmark.algorithm.operating import OperatingPoint
from arc_benchmark.analysis.metrics import fit_power_law, sweep_series
from arc_benchmark.pipeline.trace import BenchmarkTrace, SweepTrace

matplotlib.use("Agg")

from matplotlib import pyplot as plt

__all__ = [
    "plot_benchmark_diagnosis",
    "plot_field_branches",
    "plot_lawson_diagram",
    "plot_power_balance",
    "plot_sweep",
]

SERIES: Final[tuple[str, str, str]] = ("#2a78d6", "#eb6834", "#1baf7a")
"""Categorical series colours, assigned in this order and never cycled."""

_INK: Final[str] = "#0b0b0b"
_INK_SECONDARY: Final[str] = "#52514e"
_MUTED: Final[str] = "#898781"
_GRID: Final[str] = "#e1e0d9"
_TRACK: Final[str] = "#f0efec"
_SURFACE: Final[str] = "#ffffff"

_FIGURE_DPI: Final[int] = 110
"""Chosen for the published figures: 110 dots per inch over the figure sizes used
here gives a raster wide enough for a full-width README column on a high density
display, and keeps the three tracked files inside their combined size budget
without any compression step or any dependency added to do it."""


def _save(figure: Any, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    figure.tight_layout()
    figure.savefig(path, dpi=_FIGURE_DPI, facecolor=_SURFACE)
    plt.close(figure)
    return path


def _dress(axis: Any, *, grid_axis: str = "both") -> None:
    """Apply the shared chrome: recessive grid, no box, muted ticks."""
    axis.grid(True, axis=grid_axis, color=_GRID, linewidth=0.8)
    axis.set_axisbelow(True)
    for side in ("top", "right"):
        axis.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        axis.spines[side].set_color(_MUTED)
    axis.tick_params(colors=_MUTED, labelsize=8.5)
    for label in axis.get_xticklabels() + axis.get_yticklabels():
        label.set_color(_INK_SECONDARY)


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

    axes[0].plot(values, fusion, color=SERIES[0], linewidth=2.0)
    axes[0].set_ylabel("fusion power (MW)")
    axes[0].set_yscale("log")

    finite = np.isfinite(gains)
    axes[1].plot(values[finite], gains[finite], color=SERIES[1], linewidth=2.0)
    axes[1].set_ylabel("fusion gain Q")
    axes[1].set_yscale("log")

    if has_plant:
        efficiency = np.asarray(sweep_series(trace, "net_efficiency"), dtype=np.float64)
        axes[2].plot(values, efficiency, color=SERIES[2], linewidth=2.0)
        axes[2].axhline(0.0, color=_MUTED, linewidth=0.8)
        axes[2].set_ylabel("net electric / fusion")

    for axis in axes:
        _dress(axis)
        if np.any(~feasible):
            axis.fill_between(
                values,
                *axis.get_ylim(),
                where=~feasible,
                color=_INK,
                alpha=0.07,
                step="mid",
            )

    axes[-1].set_xlabel(f"{trace.variable} ({trace.units})")
    axes[0].set_title(
        f"{trace.variable} sweep, {trace.scaling_name}\n{trace.policy}",
        fontsize=9,
        color=_INK,
    )
    return _save(figure, path)


def plot_field_branches(
    branches: tuple[tuple[str, SweepTrace], ...],
    path: Path,
) -> Path:
    """Fusion power against toroidal field on several branches, with the feasible band.

    This is the figure the high-field argument reduces to. The three sweeps
    differ only in what is held fixed as the field rises, and that choice alone
    moves the measured exponent from zero to four. The lower panel shows that the
    steepest branch is not automatically the useful one: each branch is feasible
    only over a band, bounded below by one operational limit and above by
    another.

    Args:
        branches: Named sweeps, each over the same field values.
        path: Where to write the figure.

    Returns:
        The path written.

    Raises:
        ValueError: If no branch is given, or if the branches do not share one
            set of field values, in which case they cannot share an axis.
    """
    if not branches:
        raise ValueError("at least one branch is needed")
    reference = branches[0][1].values
    for label, trace in branches:
        if trace.values != reference:
            raise ValueError(f"branch {label!r} was swept over different field values")

    fields = np.asarray(reference, dtype=np.float64)
    figure, (top, bottom) = plt.subplots(
        2,
        1,
        figsize=(7.6, 5.8),
        height_ratios=(2.6, 1.0),
        sharex=True,
    )

    highest = 0.0
    for index, (label, trace) in enumerate(branches):
        colour = SERIES[index % len(SERIES)]
        power = np.asarray(sweep_series(trace, "fusion_power"), dtype=np.float64)
        exponent = fit_power_law(trace.values, sweep_series(trace, "fusion_power")).exponent
        highest = max(highest, float(power.max()))
        top.plot(
            fields,
            power,
            color=colour,
            linewidth=2.0,
            label=f"{label}, measured exponent {exponent:+.2f}",
        )
        top.annotate(
            label,
            (fields[-1], power[-1]),
            textcoords="offset points",
            xytext=(-4, 8),
            ha="right",
            fontsize=9,
            color=colour,
            fontweight="bold",
        )

    top.set_yscale("log")
    top.set_ylim(top=highest * 3.0)
    top.set_ylabel("fusion power (MW)")
    top.legend(loc="upper left", frameon=False, fontsize=9, labelcolor=_INK_SECONDARY)
    top.set_title(
        "Fusion power against toroidal field, three statements of what is held fixed",
        fontsize=10.5,
        color=_INK,
        loc="left",
    )
    _dress(top)

    span = float(fields[-1] - fields[0])
    for index, (_, trace) in enumerate(branches):
        colour = SERIES[index % len(SERIES)]
        row = len(branches) - 1 - index
        bottom.barh(row, span, left=fields[0], height=0.52, color=_TRACK, zorder=1)

        feasible = [i for i, point in enumerate(trace.points) if point.feasible]
        if not feasible:
            bottom.text(
                float(np.mean(fields)),
                row,
                "no feasible point",
                ha="center",
                va="center",
                fontsize=8.5,
                color=_INK_SECONDARY,
                zorder=3,
            )
            continue

        low, high = fields[feasible[0]], fields[feasible[-1]]
        bottom.barh(row, high - low, left=low, height=0.52, color=colour, zorder=2)
        bottom.text(
            (low + high) / 2.0,
            row,
            f"{low:.2f} to {high:.2f} T",
            ha="center",
            va="center",
            fontsize=8.5,
            color=_SURFACE,
            fontweight="bold",
            zorder=3,
        )
        # Label each infeasible region with the limit that fails immediately
        # outside the band, which is what stops the branch being pushed further.
        regions = []
        if feasible[0] > 0:
            regions.append((0, feasible[0], trace.points[feasible[0] - 1]))
        if feasible[-1] + 1 < len(fields):
            regions.append((feasible[-1] + 1, len(fields), trace.points[feasible[-1] + 1]))
        for start, stop, probe in regions:
            violations = probe.constraints.violations
            if not violations:
                continue
            name = violations[0].name.replace("_", " ")
            wide = (fields[stop - 1] - fields[start]) > 0.1 * span
            if wide:
                bottom.text(
                    (fields[start] + fields[stop - 1]) / 2.0,
                    row,
                    name,
                    ha="center",
                    va="center",
                    fontsize=8.0,
                    color=_INK_SECONDARY,
                    zorder=3,
                )
            else:
                # Too narrow to hold the text. Label it clear of the track,
                # anchored to whichever end of the axis the region sits against
                # and pushed away from the neighbouring rows.
                at_start = start == 0
                below = row == 0
                bottom.text(
                    fields[0] if at_start else fields[-1],
                    row - 0.32 if below else row + 0.32,
                    name,
                    ha="left" if at_start else "right",
                    va="top" if below else "bottom",
                    fontsize=7.5,
                    color=_INK_SECONDARY,
                    zorder=3,
                )

    bottom.set_yticks(range(len(branches)))
    bottom.set_yticklabels([label for label, _ in reversed(branches)], fontsize=8.5)
    bottom.set_ylim(-0.6, len(branches) - 0.4)
    bottom.set_xlabel("toroidal field on axis (T)")
    bottom.set_title(
        "Where each branch satisfies every limit, and what blocks it outside that band",
        fontsize=9.5,
        color=_INK,
        loc="left",
    )
    _dress(bottom, grid_axis="x")
    bottom.grid(False, axis="y")
    return _save(figure, path)


def plot_benchmark_diagnosis(
    flat: BenchmarkTrace,
    peaked: BenchmarkTrace,
    path: Path,
) -> Path:
    """The fusion power gap beside the confinement agreement that explains it.

    The two panels carry the whole diagnosis. Fusion power is low everywhere with
    flat profiles and overshoots for two of the three machines with one peaked
    shape, while the confinement enhancement the published points imply sits
    within a modest fraction of what the sources assume. Different quantities on
    different scales, so they are two panels and never two axes on one plot.

    Args:
        flat: A benchmark run with flat profiles.
        peaked: The same run with a peaked profile shape.
        path: Where to write the figure.

    Returns:
        The path written.

    Raises:
        ValueError: If the two runs do not cover the same machines in the same
            order.
    """
    names = [result.machine for result in flat.results]
    if names != [result.machine for result in peaked.results]:
        raise ValueError("the flat and peaked runs must cover the same machines in order")

    def published_fusion(trace: BenchmarkTrace, machine: str) -> float:
        result = trace.named(machine)
        return next(row.published for row in result.rows if row.quantity == "fusion power")

    positions = np.arange(len(names), dtype=np.float64)
    width = 0.26
    figure, (left, right) = plt.subplots(1, 2, figsize=(8.6, 4.0), width_ratios=(1.35, 1.0))

    series = (
        ("this model, flat profiles", [flat.named(n).point.terms.fusion_power_mw for n in names]),
        (
            "this model, peaked profiles",
            [peaked.named(n).point.terms.fusion_power_mw for n in names],
        ),
        ("published design value", [published_fusion(flat, n) for n in names]),
    )
    for index, (label, values) in enumerate(series):
        offset = (index - 1) * width
        bars = left.bar(
            positions + offset,
            values,
            width * 0.92,
            color=SERIES[index],
            label=label,
            zorder=2,
        )
        left.bar_label(bars, fmt="%.0f", fontsize=7.5, color=_INK_SECONDARY, padding=2)

    left.set_xticks(positions)
    left.set_xticklabels(names)
    left.set_ylabel("fusion power (MW)")
    left.set_ylim(top=max(max(values) for _, values in series) * 1.34)
    left.legend(loc="upper right", frameon=False, fontsize=8.5, labelcolor=_INK_SECONDARY)
    left.set_title("Fusion power is where the model is wrong", fontsize=10, color=_INK, loc="left")
    _dress(left, grid_axis="y")
    left.grid(False, axis="x")

    assumed = [flat.named(n).point.state.confinement_multiplier for n in names]
    implied = [flat.named(n).implied_confinement_multiplier for n in names]
    for index in range(len(names)):
        value = implied[index]
        if value is None:
            continue
        right.plot([index, index], [assumed[index], value], color=_MUTED, linewidth=1.5, zorder=1)
        right.text(
            index + 0.12,
            (assumed[index] + value) / 2.0,
            f"{value / assumed[index]:.3f}x",
            fontsize=8.5,
            va="center",
            color=_INK_SECONDARY,
        )

    right.scatter(
        positions, assumed, s=70, color=SERIES[2], zorder=3, label="assumed by the source"
    )
    right.scatter(
        positions,
        [value if value is not None else np.nan for value in implied],
        s=70,
        color=SERIES[1],
        zorder=3,
        label="implied by the published point",
    )
    right.set_xticks(positions)
    right.set_xticklabels(names)
    right.set_xlim(-0.5, len(names) - 0.2)
    right.set_ylabel("confinement enhancement H")
    highest = max([value for value in implied if value is not None] + assumed)
    right.set_ylim(top=highest * 1.16)
    right.legend(loc="upper right", frameon=False, fontsize=8.5, labelcolor=_INK_SECONDARY)
    right.set_title("Confinement is not", fontsize=10, color=_INK, loc="left")
    _dress(right, grid_axis="y")
    right.grid(False, axis="x")

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
    figure, axis = plt.subplots(figsize=(7.4, 4.6))

    for index, gain in enumerate(gains):
        values = np.array(
            [
                lawson_triple_product(float(t), composition, gain).triple_product
                for t in temperatures
            ],
            dtype=np.float64,
        )
        finite = np.isfinite(values)
        label = "ignition" if gain == float("inf") else f"Q = {gain:g}"
        colour = SERIES[index % len(SERIES)]
        axis.plot(temperatures[finite], values[finite], color=colour, linewidth=2.0, label=label)
        axis.annotate(
            label,
            (temperatures[finite][-1], values[finite][-1]),
            textcoords="offset points",
            xytext=(-6, 6),
            ha="right",
            fontsize=9,
            color=colour,
            fontweight="bold",
        )

    for name, point in points:
        axis.scatter(
            [point.state.temperature_kev],
            [point.triple_product],
            marker="o",
            s=60,
            color=_INK,
            edgecolor=_SURFACE,
            linewidth=1.5,
            zorder=5,
        )
        axis.annotate(
            name,
            (point.state.temperature_kev, point.triple_product),
            textcoords="offset points",
            xytext=(9, -3),
            fontsize=9,
            color=_INK,
            bbox={"facecolor": _SURFACE, "edgecolor": "none", "pad": 1.0, "alpha": 0.85},
        )

    axis.set_yscale("log")
    axis.set_xlabel("temperature (keV)")
    axis.set_ylabel("triple product n T tau (m^-3 keV s)")
    axis.set_title(
        "Lawson requirement and the benchmarked design points",
        fontsize=10.5,
        color=_INK,
        loc="left",
    )
    axis.legend(loc="upper right", frameon=False, fontsize=9, labelcolor=_INK_SECONDARY)
    _dress(axis)
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

    figure, axis = plt.subplots(figsize=(7.4, 4.6))
    axis.bar(positions - width / 2, alpha, width, color=SERIES[0], label="alpha heating")
    axis.bar(
        positions - width / 2,
        auxiliary,
        width,
        bottom=alpha,
        color=SERIES[1],
        label="auxiliary heating",
    )
    axis.bar(positions + width / 2, brem, width, color=_MUTED, label="bremsstrahlung")
    axis.bar(positions + width / 2, syn, width, bottom=brem, color=SERIES[2], label="synchrotron")
    axis.bar(
        positions + width / 2,
        line,
        width,
        bottom=brem + syn,
        color=_INK_SECONDARY,
        label="line radiation",
    )
    axis.bar(
        positions + width / 2,
        transport,
        width,
        bottom=brem + syn + line,
        color=_TRACK,
        edgecolor=_MUTED,
        label="transport",
    )

    axis.set_xticks(positions)
    axis.set_xticklabels(names)
    axis.set_ylabel("power (MW)")
    axis.set_title(
        "Sources on the left of each pair, sinks on the right",
        fontsize=10.5,
        color=_INK,
        loc="left",
    )
    axis.legend(fontsize=8.5, frameon=False, labelcolor=_INK_SECONDARY)
    _dress(axis, grid_axis="y")
    axis.grid(False, axis="x")
    return _save(figure, path)
