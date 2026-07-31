"""Net electrical efficiency: where the recirculating power goes and what recovers it.

This is the study the repository exists for. The plasma gain is an input to it,
not the answer: a plant with a high plasma gain and a poor current drive
efficiency can still consume more electricity than it produces.
"""

from __future__ import annotations

import argparse
import dataclasses
from collections.abc import Sequence
from pathlib import Path

import numpy as np

from arc_benchmark.algorithm.operating import solve_operating_point
from arc_benchmark.analysis.figures import plot_sweep
from arc_benchmark.analysis.metrics import sweep_series
from arc_benchmark.analysis.report import sweep_lines
from arc_benchmark.analysis.tables import sweep_frame
from arc_benchmark.model.confinement import IPB98Y2
from arc_benchmark.model.plant import evaluate_plant
from arc_benchmark.pipeline.machines import machine
from arc_benchmark.pipeline.sweep import density_sweep


def _parse(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--points", type=int, default=41, help="points in the density sweep")
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
    parameters = case.plant
    assert parameters is not None, "the ARC case defines a plant"

    point = solve_operating_point(base, IPB98Y2)
    plant = evaluate_plant(
        fusion_power_mw=point.terms.fusion_power_mw,
        alpha_power_mw=point.terms.alpha_power_mw,
        neutron_power_mw=point.terms.neutron_power_mw,
        auxiliary_power_mw=max(point.auxiliary_power_mw, 0.0),
        parameters=parameters,
    )

    print("Electrical accounting at the base design point")
    print(f"  fusion power                 {point.terms.fusion_power_mw:10.2f} MW")
    print(f"  neutron power                {point.terms.neutron_power_mw:10.2f} MW")
    print(f"  blanket multiplication       {parameters.blanket_multiplication:10.2f}")
    print(f"  thermal power                {plant.thermal_power_mw:10.2f} MW")
    print(f"  thermal efficiency           {parameters.thermal_efficiency:10.2f}")
    print(f"  gross electric               {plant.gross_electric_mw:10.2f} MWe")
    print("  recirculating draws")
    print(f"    heating and current drive  {plant.heating_draw_mw:10.2f} MWe")
    print(f"    cryoplant                  {plant.cryoplant_mw:10.2f} MWe")
    print(f"    balance of plant           {plant.balance_of_plant_mw:10.2f} MWe")
    print(f"    tritium and site           {plant.tritium_and_auxiliary_mw:10.2f} MWe")
    print(f"    total                      {plant.recirculating_mw:10.2f} MWe")
    print(f"  net electric                 {plant.net_electric_mw:10.2f} MWe")
    print(f"  recirculating fraction       {plant.recirculating_fraction:10.3f}")
    print(f"  net efficiency               {plant.net_efficiency:10.4f}")
    print(f"  engineering gain             {plant.engineering_gain:10.3f}")

    densities = np.linspace(0.6e20, 2.4e20, points)
    trace = density_sweep(base, densities, IPB98Y2, plant=parameters)
    print()
    for line in sweep_lines(trace, stride=max(1, points // 9)):
        print(line)

    efficiencies = np.asarray(sweep_series(trace, "net_efficiency"))
    best = int(np.argmax(efficiencies))
    print(
        f"best net efficiency {efficiencies[best]:.4f} at density "
        f"{trace.values[best]:.4e} m^-3, feasible {trace.points[best].feasible}, "
        f"binding {trace.points[best].binding_constraint}"
    )
    feasible = trace.feasible_points
    if feasible:
        best_feasible = max(
            feasible, key=lambda p: p.plant.net_efficiency if p.plant is not None else -1.0
        )
        assert best_feasible.plant is not None
        print(
            f"best feasible net efficiency {best_feasible.plant.net_efficiency:.4f} at density "
            f"{best_feasible.value:.4e} m^-3"
        )

    print()
    print("The five feasible densities with the highest net efficiency")
    frame = sweep_frame(trace)
    ranked = frame[frame["feasible"]].sort_values("net_efficiency", ascending=False)
    print(
        ranked[
            [
                "electron density",
                "fusion_power_mw",
                "auxiliary_power_mw",
                "net_electric_mw",
                "net_efficiency",
                "binding_constraint",
            ]
        ]
        .head(5)
        .to_string(index=False)
    )

    print()
    print("Sensitivity of net electric power to one engineering parameter at a time")
    print(f"  {'parameter':34} {'value':>8} {'net MWe':>10} {'eta_net':>9}")
    variations = (
        ("thermal_efficiency", (0.33, 0.40, 0.45, 0.50)),
        ("heating_wallplug_efficiency", (0.35, 0.45, 0.55, 0.70)),
        ("blanket_multiplication", (1.10, 1.20, 1.30, 1.40)),
        ("cryoplant_mw", (2.0, 10.0, 25.0, 50.0)),
    )
    for name, values in variations:
        for value in values:
            trial = dataclasses.replace(parameters, **{name: value})
            result = evaluate_plant(
                fusion_power_mw=point.terms.fusion_power_mw,
                alpha_power_mw=point.terms.alpha_power_mw,
                neutron_power_mw=point.terms.neutron_power_mw,
                auxiliary_power_mw=max(point.auxiliary_power_mw, 0.0),
                parameters=trial,
            )
            print(
                f"  {name:34} {value:8.2f} {result.net_electric_mw:10.2f} "
                f"{result.net_efficiency:9.4f}"
            )

    print()
    print("Sensitivity of the plasma solution to the first wall reflectivity")
    print(f"  {'reflectivity':>12} {'P_syn MW':>10} {'P_aux MW':>10} {'Q':>8} {'net MWe':>10}")
    for reflectivity in (0.6, 0.8, 0.9, 0.95, 0.98):
        trial_state = dataclasses.replace(base, wall_reflectivity=reflectivity)
        trial_point = solve_operating_point(trial_state, IPB98Y2)
        trial_plant = evaluate_plant(
            fusion_power_mw=trial_point.terms.fusion_power_mw,
            alpha_power_mw=trial_point.terms.alpha_power_mw,
            neutron_power_mw=trial_point.terms.neutron_power_mw,
            auxiliary_power_mw=max(trial_point.auxiliary_power_mw, 0.0),
            parameters=parameters,
        )
        print(
            f"  {reflectivity:12.2f} {trial_point.terms.synchrotron_mw:10.2f} "
            f"{trial_point.auxiliary_power_mw:10.2f} {trial_point.fusion_gain:8.2f} "
            f"{trial_plant.net_electric_mw:10.2f}"
        )

    if not args.quick:
        path = plot_sweep(trace, args.figure_dir / "efficiency_density_sweep.png")
        print()
        print(f"wrote {path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
