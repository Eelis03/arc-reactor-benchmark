"""Measure how fusion power, gain, and net efficiency scale with toroidal field.

The high-field pathway is an argument about an exponent, so this script measures
the exponent under three different statements of what is held fixed, and reports
which operational limit binds first in each case.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from arc_benchmark.analysis.figures import plot_sweep
from arc_benchmark.analysis.metrics import binding_constraint_counts, fit_power_law, sweep_series
from arc_benchmark.analysis.report import sweep_lines
from arc_benchmark.model.confinement import IPB98Y2
from arc_benchmark.pipeline.machines import machine
from arc_benchmark.pipeline.sweep import SweepInvariant, field_sweep, radius_sweep


def _parse(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field-min", type=float, default=4.0, help="lowest field in tesla")
    parser.add_argument("--field-max", type=float, default=14.0, help="highest field in tesla")
    parser.add_argument("--points", type=int, default=41, help="points per sweep")
    parser.add_argument("--radius-min", type=float, default=2.0, help="lowest major radius")
    parser.add_argument("--radius-max", type=float, default=8.0, help="highest major radius")
    parser.add_argument(
        "--figure-dir", type=Path, default=Path("figures"), help="figure output directory"
    )
    parser.add_argument("--quick", action="store_true", help="fewer points and no figures")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the example. Returns a process exit code."""
    args = _parse(argv)
    points = 9 if args.quick else args.points
    case = machine("ARC")
    base = case.state
    plant = case.plant

    fields = np.linspace(args.field_min, args.field_max, points)
    stride = max(1, points // 9)

    expected = {
        SweepInvariant.FIXED_DENSITY: (
            "0, because fusion power at fixed density, temperature, and volume does "
            "not depend on the field at all"
        ),
        SweepInvariant.FIXED_GREENWALD_FRACTION: (
            "2, because the density tracks the current, the current tracks the field "
            "at fixed safety factor, and fusion power goes as the square of density"
        ),
        SweepInvariant.FIXED_BETA: (
            "4, because holding beta at fixed temperature makes the density go as the "
            "square of the field, and fusion power as the square of that"
        ),
    }

    for invariant in SweepInvariant:
        trace = field_sweep(base, fields, invariant, IPB98Y2, plant=plant)
        print()
        for line in sweep_lines(trace, stride=stride):
            print(line)

        fusion_fit = fit_power_law(trace.values, sweep_series(trace, "fusion_power"))
        print(
            f"fusion power exponent in field: measured {fusion_fit.exponent:+.4f}, "
            f"r squared {fusion_fit.r_squared:.6f}"
        )
        print(f"  analytic expectation: {expected[invariant]}")

        gains = sweep_series(trace, "gain")
        if all(np.isfinite(gains)):
            gain_fit = fit_power_law(trace.values, gains)
            print(
                f"fusion gain exponent in field: measured {gain_fit.exponent:+.4f}, "
                f"r squared {gain_fit.r_squared:.6f}"
            )

        tau_fit = fit_power_law(trace.values, sweep_series(trace, "confinement_time"))
        print(f"confinement time exponent in field: measured {tau_fit.exponent:+.4f}")

        if not args.quick:
            name = invariant.name.lower()
            path = plot_sweep(trace, args.figure_dir / f"field_sweep_{name}.png")
            print(f"wrote {path}")

    radii = np.linspace(args.radius_min, args.radius_max, points)
    radius_trace = radius_sweep(
        base, radii, SweepInvariant.FIXED_GREENWALD_FRACTION, IPB98Y2, plant=plant
    )
    print()
    for line in sweep_lines(radius_trace, stride=stride):
        print(line)
    radius_fit = fit_power_law(radius_trace.values, sweep_series(radius_trace, "fusion_power"))
    print(
        f"fusion power exponent in major radius: measured {radius_fit.exponent:+.4f}, "
        f"r squared {radius_fit.r_squared:.6f}"
    )
    print(
        "  analytic expectation: 1, because volume goes as the cube of the radius while "
        "the Greenwald density falls as the inverse radius at fixed safety factor"
    )
    print(f"binding constraints: {binding_constraint_counts(radius_trace)}")

    if not args.quick:
        path = plot_sweep(radius_trace, args.figure_dir / "radius_sweep.png")
        print(f"wrote {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
