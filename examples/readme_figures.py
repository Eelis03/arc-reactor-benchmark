"""Regenerate the three figures the README publishes, and print what is in them.

This is the one command that produces everything under ``docs/figures``. The
figures are snapshots: they are written by this script, committed as files, and
never compared byte for byte by continuous integration, because matplotlib output
is not byte reproducible across platforms and font stacks.

Every number printed here is the number the corresponding figure draws, so a
caption can be checked against this output rather than against the picture.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from arc_benchmark.algorithm.balance import PlasmaComposition
from arc_benchmark.analysis.figures import (
    plot_benchmark_diagnosis,
    plot_field_branches,
    plot_lawson_diagram,
)
from arc_benchmark.analysis.metrics import binding_constraint_counts, fit_power_law, sweep_series
from arc_benchmark.model.confinement import IPB98Y2
from arc_benchmark.model.profiles import ProfileShape
from arc_benchmark.pipeline.benchmark import run_benchmark
from arc_benchmark.pipeline.machines import machine
from arc_benchmark.pipeline.sweep import SweepInvariant, field_sweep

_FIGURE_BUDGET_BYTES = 250 * 1024
"""The three tracked figures have to fit in this together. Checked and reported
here rather than discovered when a repository grows a megabyte of images."""

_BRANCH_LABELS = {
    SweepInvariant.FIXED_DENSITY: "fixed electron density",
    SweepInvariant.FIXED_GREENWALD_FRACTION: "fixed Greenwald fraction",
    SweepInvariant.FIXED_BETA: "fixed toroidal beta",
}


def _parse(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=41, help="points per field sweep")
    parser.add_argument("--field-min", type=float, default=4.0, help="lowest field in tesla")
    parser.add_argument("--field-max", type=float, default=14.0, help="highest field in tesla")
    parser.add_argument(
        "--density-exponent", type=float, default=0.4, help="peaked case density exponent"
    )
    parser.add_argument(
        "--temperature-exponent", type=float, default=1.0, help="peaked case temperature exponent"
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("docs") / "figures",
        help="directory the tracked figures are written to",
    )
    parser.add_argument("--quick", action="store_true", help="fewer points and no figures")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the example. Returns a process exit code."""
    args = _parse(argv)
    points = 9 if args.quick else args.points
    case = machine("ARC")
    fields = np.linspace(args.field_min, args.field_max, points)

    branches = tuple(
        (_BRANCH_LABELS[invariant], field_sweep(case.state, fields, invariant, IPB98Y2))
        for invariant in SweepInvariant
    )

    print("Figure 1, fusion power against toroidal field on three branches")
    for label, trace in branches:
        fit = fit_power_law(trace.values, sweep_series(trace, "fusion_power"))
        feasible = trace.feasible_points
        band = (
            f"{feasible[0].value:.2f} to {feasible[-1].value:.2f} T"
            if feasible
            else "no feasible point"
        )
        print(
            f"  {label:26} exponent {fit.exponent:+.4f}, r squared {fit.r_squared:.6f}, "
            f"feasible {len(feasible)} of {len(trace.points)}, {band}"
        )
        print(f"    binding constraint tally: {binding_constraint_counts(trace)}")

    flat = run_benchmark(IPB98Y2)
    peaked = run_benchmark(
        IPB98Y2, profile=ProfileShape(args.density_exponent, args.temperature_exponent)
    )

    print()
    print("Figure 2, the Lawson requirement and the three benchmarked points")
    for result in flat.results:
        point = result.point
        print(
            f"  {result.machine:6} triple product {point.triple_product:10.4e} m^-3 keV s "
            f"at {point.state.temperature_kev:5.2f} keV"
        )

    print()
    print("Figure 3, the fusion power gap and the implied confinement enhancement")
    print(f"  {'machine':8} {'flat MW':>9} {'peaked MW':>10} {'published MW':>13} "
          f"{'H assumed':>10} {'H implied':>10} {'ratio':>7}")
    for result in flat.results:
        name = result.machine
        published = next(row.published for row in result.rows if row.quantity == "fusion power")
        assumed = result.point.state.confinement_multiplier
        implied = result.implied_confinement_multiplier
        ratio = implied / assumed if implied is not None else float("nan")
        print(
            f"  {name:8} {result.point.terms.fusion_power_mw:9.1f} "
            f"{peaked.named(name).point.terms.fusion_power_mw:10.1f} {published:13.1f} "
            f"{assumed:10.3f} {implied if implied is not None else float('nan'):10.3f} "
            f"{ratio:7.3f}"
        )

    if args.quick:
        return 0

    written = (
        plot_field_branches(branches, args.figure_dir / "field_scaling.png"),
        plot_lawson_diagram(
            PlasmaComposition(helium_fraction=0.05),
            tuple((result.machine, result.point) for result in flat.results),
            args.figure_dir / "lawson.png",
        ),
        plot_benchmark_diagnosis(flat, peaked, args.figure_dir / "benchmark_diagnosis.png"),
    )

    print()
    total = 0
    for path in written:
        size = path.stat().st_size
        total += size
        print(f"wrote {path} ({size} bytes)")
    print(f"total {total} bytes against a budget of {_FIGURE_BUDGET_BYTES}")
    return 0 if total <= _FIGURE_BUDGET_BYTES else 1


if __name__ == "__main__":
    raise SystemExit(main())
