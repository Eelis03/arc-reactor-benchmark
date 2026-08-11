"""How much of the answer is the confinement scaling rather than the physics.

Runs the same benchmark and the same field sweep under three confinement
scalings, and measures the amplification factor by which an error in the
confinement time becomes an error in the required heating power.
"""

from __future__ import annotations

import argparse
import dataclasses
from collections.abc import Sequence

import numpy as np

from arc_benchmark.algorithm.balance import LossPowerConvention
from arc_benchmark.algorithm.constraints import ConstraintLimits
from arc_benchmark.algorithm.operating import solve_operating_point
from arc_benchmark.analysis.metrics import fit_power_law, sweep_series
from arc_benchmark.model.confinement import CONFINEMENT_SCALINGS
from arc_benchmark.pipeline.benchmark import run_benchmark
from arc_benchmark.pipeline.machines import machine
from arc_benchmark.pipeline.sweep import SweepInvariant, field_sweep


def _parse(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=41, help="points in the field sweep")
    parser.add_argument("--quick", action="store_true", help="fewer points")
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the example. Returns a process exit code."""
    args = _parse(argv)
    points = 9 if args.quick else args.points
    case = machine("ARC")

    print("Power degradation exponent and its amplification into the loss power")
    print(f"  {'scaling':12} {'mode':>5} {'alpha_P':>9} {'1/(1-alpha_P)':>15}")
    for name, scaling in CONFINEMENT_SCALINGS.items():
        alpha = scaling.power_degradation
        print(f"  {name:12} {scaling.confinement_mode:>5} {alpha:9.3f} {1.0 / (1.0 - alpha):15.3f}")
    print(
        "  The loss power at a solved steady state goes as the stored energy to the "
        "power 1 / (1 - alpha_P), so an error in the confinement time is amplified by "
        "that exponent before it reaches the required heating power."
    )

    print()
    print("The same three design points under each scaling")
    for name, scaling in CONFINEMENT_SCALINGS.items():
        limits = ConstraintLimits(require_h_mode=scaling.confinement_mode == "H")
        trace = run_benchmark(scaling, limits=limits)
        print(f"  {name} ({scaling.confinement_mode}-mode fit)")
        print(
            f"    {'machine':8} {'tau_E s':>9} {'P_aux MW':>10} {'Q':>8} "
            f"{'implied H':>10} {'binding':>14}"
        )
        for result in trace.results:
            implied = result.implied_confinement_multiplier
            implied_text = "n/a" if implied is None else f"{implied:.3f}"
            print(
                f"    {result.machine:8} {result.point.confinement_time_s:9.3f} "
                f"{result.point.auxiliary_power_mw:10.2f} {result.point.fusion_gain:8.2f} "
                f"{implied_text:>10} {result.constraints.binding.name:>14}"
            )

    print()
    print("Direct sensitivity of the ARC point to the assumed H factor")
    print(f"  {'H':>6} {'tau_E s':>9} {'P_loss MW':>11} {'P_aux MW':>10} {'Q':>8}")
    reference = solve_operating_point(case.state, CONFINEMENT_SCALINGS["IPB98(y,2)"])
    for multiplier in (1.2, 1.5, 1.8, 2.1, 2.5, 3.0):
        trial = dataclasses.replace(case.state, confinement_multiplier=multiplier)
        solved = solve_operating_point(trial, CONFINEMENT_SCALINGS["IPB98(y,2)"])
        print(
            f"  {multiplier:6.2f} {solved.confinement_time_s:9.3f} "
            f"{solved.loss_power_mw:11.2f} {solved.auxiliary_power_mw:10.2f} "
            f"{solved.fusion_gain:8.2f}"
        )
    print(
        f"  reference point at H = {case.state.confinement_multiplier}: "
        f"Q = {reference.fusion_gain:.2f}"
    )

    print()
    print("Loss power convention, which changes whether radiation is subtracted")
    print(f"  {'convention':12} {'P_loss MW':>11} {'P_aux MW':>10} {'Q':>8}")
    for convention in LossPowerConvention:
        solved = solve_operating_point(case.state, CONFINEMENT_SCALINGS["IPB98(y,2)"], convention)
        print(
            f"  {convention.value:12} {solved.loss_power_mw:11.2f} "
            f"{solved.auxiliary_power_mw:10.2f} {solved.fusion_gain:8.2f}"
        )

    print()
    print("Field scaling exponent of fusion power under each scaling, at fixed beta")
    fields = np.linspace(4.0, 14.0, points)
    for name, scaling in CONFINEMENT_SCALINGS.items():
        limits = ConstraintLimits(require_h_mode=scaling.confinement_mode == "H")
        trace = field_sweep(
            case.state,
            fields,
            SweepInvariant.FIXED_BETA,
            scaling,
            limits=limits,
            plant=case.plant,
        )
        fusion_fit = fit_power_law(trace.values, sweep_series(trace, "fusion_power"))
        gains = sweep_series(trace, "gain")
        gain_text = "not finite everywhere"
        if all(np.isfinite(gains)):
            gain_text = f"{fit_power_law(trace.values, gains).exponent:+.4f}"
        feasible = len(trace.feasible_points)
        print(
            f"  {name:12} fusion {fusion_fit.exponent:+.4f}  gain {gain_text:>22}  "
            f"feasible {feasible}/{len(trace.points)}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
