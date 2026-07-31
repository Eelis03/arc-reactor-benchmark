"""Solve and print the zero-dimensional power balance of the ARC design point.

Wiring only. Every number printed here is computed in the library.
"""

from __future__ import annotations

import argparse
from collections.abc import Sequence
from pathlib import Path

from arc_benchmark.algorithm.constraints import evaluate_constraints
from arc_benchmark.algorithm.lawson import lawson_triple_product, optimum_lawson_temperature
from arc_benchmark.algorithm.operating import solve_ignition_temperature, solve_operating_point
from arc_benchmark.analysis.figures import plot_lawson_diagram, plot_power_balance
from arc_benchmark.analysis.report import constraint_lines, power_balance_lines
from arc_benchmark.model.confinement import IPB98Y2
from arc_benchmark.pipeline.machines import MACHINES, machine


def _parse(argv: Sequence[str] | None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--machine",
        default="ARC",
        choices=sorted(MACHINES),
        help="which published design point to solve",
    )
    parser.add_argument(
        "--figure-dir",
        type=Path,
        default=Path("figures"),
        help="directory the figures are written to",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="skip the figures, for a short run",
    )
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the example. Returns a process exit code."""
    args = _parse(argv)
    case = machine(args.machine)

    point = solve_operating_point(case.state, IPB98Y2)
    report = evaluate_constraints(point)

    print(f"{case.name}: {case.description}")
    print(f"source: {case.source}")
    print()
    for line in power_balance_lines(point):
        print(line)
    print()
    for line in constraint_lines(report):
        print(line)

    print()
    composition = case.state.composition
    required = lawson_triple_product(case.state.temperature_kev, composition)
    optimum = optimum_lawson_temperature(composition)
    achieved = point.triple_product
    print("Lawson condition, evaluated for this composition")
    print(f"  achieved triple product        {achieved:10.4e} m^-3 keV s")
    print(
        f"  ignition requirement at {case.state.temperature_kev:5.1f} keV "
        f"{required.triple_product:10.4e} m^-3 keV s"
    )
    print(f"  ratio achieved to required     {achieved / required.triple_product:10.3f}")
    print(
        f"  minimum requirement            {optimum.triple_product:10.4e} m^-3 keV s "
        f"at {optimum.temperature_kev:.2f} keV"
    )

    ignition_temperature = solve_ignition_temperature(case.state, IPB98Y2)
    if ignition_temperature is None:
        print("  this state does not ignite anywhere in the searched temperature bracket")
    else:
        print(f"  ignites at                     {ignition_temperature:10.2f} keV")

    if not args.quick:
        balance_path = plot_power_balance(((case.name, point),), args.figure_dir / "balance.png")
        lawson_path = plot_lawson_diagram(
            composition, ((case.name, point),), args.figure_dir / "lawson.png"
        )
        print()
        print(f"wrote {balance_path}")
        print(f"wrote {lawson_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
