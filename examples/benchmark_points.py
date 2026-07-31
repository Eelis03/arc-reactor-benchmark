"""Reproduce published design points and report every discrepancy.

Runs the benchmark twice: once with flat profiles, which is what a genuinely
zero-dimensional model produces, and once with a parabolic profile shape, so that
the size of the profile effect can be read off rather than argued about.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from arc_benchmark.algorithm.balance import PlasmaComposition
from arc_benchmark.analysis.figures import plot_lawson_diagram, plot_power_balance
from arc_benchmark.analysis.report import benchmark_lines
from arc_benchmark.analysis.tables import benchmark_frame
from arc_benchmark.model.confinement import IPB98Y2
from arc_benchmark.model.profiles import ProfileShape
from arc_benchmark.pipeline.benchmark import run_benchmark


def _parse(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--density-exponent",
        type=float,
        default=0.4,
        help="parabolic density profile exponent for the peaked case",
    )
    parser.add_argument(
        "--temperature-exponent",
        type=float,
        default=1.0,
        help="parabolic temperature profile exponent for the peaked case",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("figures"),
        help="directory the figures are written to",
    )
    parser.add_argument("--quick", action="store_true", help="skip the figures")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the example. Returns a process exit code."""
    args = _parse(argv)
    peaked = ProfileShape(args.density_exponent, args.temperature_exponent)

    flat_trace = run_benchmark(IPB98Y2)
    print("FLAT PROFILES, the zero-dimensional answer")
    for line in benchmark_lines(flat_trace):
        print(line)

    print()
    print(
        f"PARABOLIC PROFILES, density exponent {args.density_exponent}, "
        f"temperature exponent {args.temperature_exponent}"
    )
    print(
        f"  fusion enhancement at 14 keV {peaked.fusion_factor(14.0):.3f}, "
        f"stored energy {peaked.stored_energy_factor():.3f}, "
        f"bremsstrahlung {peaked.bremsstrahlung_factor():.3f}"
    )
    peaked_trace = run_benchmark(IPB98Y2, profile=peaked)
    for line in benchmark_lines(peaked_trace):
        print(line)

    print()
    print("Fusion power, flat against peaked, against published")
    print(f"  {'machine':8} {'flat MW':>10} {'peaked MW':>10} {'published MW':>13}")
    for flat, peak in zip(flat_trace.results, peaked_trace.results, strict=True):
        published = next(
            (row.published for row in flat.rows if row.quantity == "fusion power"), float("nan")
        )
        print(
            f"  {flat.machine:8} {flat.point.terms.fusion_power_mw:10.1f} "
            f"{peak.point.terms.fusion_power_mw:10.1f} {published:13.1f}"
        )

    print()
    print("The same flat-profile comparison as a data frame, sorted by discrepancy")
    frame = benchmark_frame(flat_trace)
    frame = frame.assign(absolute_relative_error=frame["relative_error"].abs())
    ranked = frame.sort_values("absolute_relative_error", ascending=False)
    print(
        ranked[["machine", "quantity", "computed", "published", "relative_error"]]
        .head(10)
        .to_string(index=False)
    )

    if not args.quick:
        labelled = tuple((result.machine, result.point) for result in flat_trace.results)
        balance_path = plot_power_balance(labelled, args.figure_dir / "benchmark_balance.png")
        lawson_path = plot_lawson_diagram(
            PlasmaComposition(helium_fraction=0.05),
            labelled,
            args.figure_dir / "benchmark_lawson.png",
        )
        print()
        print(f"wrote {balance_path}")
        print(f"wrote {lawson_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
