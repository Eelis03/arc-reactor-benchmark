"""Tier two: a recorded reference run, pinned with tolerances derived from the maths.

What is pinned here, and what deliberately is not.

Every operating point in this package is the result of a closed-form inversion.
The steady-state balance reduces to ``P_loss = (W / tau_1)**(1 / (1 - alpha))``,
which is a handful of arithmetic operations with no iteration, no tolerance, and
no convergence test, so the same inputs produce the same bits on any machine that
implements IEEE 754. Those values are safe to pin tightly and are pinned tightly.

The two iterative solves in the package are Brent's method for an ignition
temperature and a bounded Brent minimisation for the Lawson optimum. Both check
convergence before returning, and both are pinned only to the tolerance they were
asked to converge to, not to the digits they happen to produce. An unconverged
result is never pinned: :func:`solve_ignition_temperature` returns ``None``
rather than a number when the bracket contains no sign change, and that ``None``
is what the tests assert, because a state that does not ignite is a result to
report and not a number to record.

Sweep-derived quantities are quantised by the sweep resolution. The field at
which a constraint first binds cannot be resolved more finely than one sweep
step, so it is pinned to one step and not to the floating-point value that step
happens to land on. Counts are integers and are pinned exactly.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from arc_benchmark.algorithm.balance import PlasmaComposition
from arc_benchmark.algorithm.lawson import optimum_lawson_temperature
from arc_benchmark.algorithm.operating import solve_operating_point
from arc_benchmark.analysis.metrics import binding_constraint_counts, fit_power_law, sweep_series
from arc_benchmark.model.confinement import IPB98Y2, ITER89P, PETTY08
from arc_benchmark.pipeline.benchmark import run_benchmark
from arc_benchmark.pipeline.machines import machine
from arc_benchmark.pipeline.sweep import SweepInvariant, density_sweep, field_sweep, radius_sweep

# Tolerance for a closed-form quantity. Each of these values is reached by fewer
# than a hundred floating-point operations, so the accumulated relative error is
# bounded by roughly 100 * 2.2e-16, which is 2e-14. The pin is set at 1e-10,
# four orders above that bound. It is derived from the operation count and the
# machine epsilon, not from any difference observed between two machines.
_CLOSED_FORM_RTOL = 1.0e-10

_OPERATING_POINTS: dict[str, dict[str, float]] = {
    "ARC": {
        "fusion_power_mw": 355.10448888932456,
        "alpha_power_mw": 71.0237241909936,
        "auxiliary_power_mw": 112.2237734788986,
        "fusion_gain": 3.164253686016847,
        "confinement_time_s": 1.0753236982306174,
        "loss_power_mw": 121.39109418739574,
        "stored_energy_mj": 130.5347203338516,
        "bremsstrahlung_mw": 6.061181006896073,
        "synchrotron_mw": 55.79522247560037,
        "line_radiation_mw": 0.0,
        "triple_product": 1.9570891307797237e21,
        "reactivity_m3_s": 2.405841927609959e-22,
        "implied_h": 2.0347116850131206,
    },
    "ITER": {
        "fusion_power_mw": 335.6296816034005,
        "alpha_power_mw": 67.12860772633152,
        "auxiliary_power_mw": 91.68566208176225,
        "fusion_gain": 3.660656137315091,
        "confinement_time_s": 3.0901398857177447,
        "loss_power_mw": 108.1740866233505,
        "stored_energy_mj": 334.2730596759017,
        "bremsstrahlung_mw": 18.87445877444364,
        "synchrotron_mw": 25.10808406550079,
        "line_radiation_mw": 6.657640344798839,
        "triple_product": 2.7193230994316156e21,
        "reactivity_m3_s": 8.127145026735505e-23,
        "implied_h": 1.0266855949666702,
    },
    "SPARC": {
        "fusion_power_mw": 66.41901735091318,
        "alpha_power_mw": 13.284332124673101,
        "auxiliary_power_mw": 36.185306953197355,
        "fusion_gain": 1.8355244971894416,
        "confinement_time_s": 0.787792534559029,
        "loss_power_mw": 31.787617530027354,
        "stored_energy_mj": 25.04204778157327,
        "bremsstrahlung_mw": 3.5556950642507466,
        "synchrotron_mw": 14.126326483592356,
        "line_radiation_mw": 0.0,
        "triple_product": 1.7827745057070826e21,
        "reactivity_m3_s": 4.750036924137603e-23,
        "implied_h": 1.1318312233685428,
    },
}

_CONSTRAINT_UTILISATION: dict[str, dict[str, float]] = {
    "ARC": {
        "greenwald": 0.6685832765614675,
        "troyon": 0.803700493704306,
        "safety_factor": 0.28139027255823734,
        "bootstrap": 0.6468695732488934,
        "lh_threshold": 0.44883710865722587,
    },
    "ITER": {
        "greenwald": 0.8377580409572781,
        "troyon": 0.6046784203343877,
        "safety_factor": 0.6674813134460509,
        "bootstrap": 0.22410809439577942,
        "lh_threshold": 0.7878754134384138,
    },
    "SPARC": {
        "greenwald": 0.36369893158265565,
        "troyon": 0.3443013890375086,
        "safety_factor": 0.42802015690198025,
        "bootstrap": 0.17404772782069136,
        "lh_threshold": 1.279672863976897,
    },
}

# Verdicts are categorical and are pinned exactly. SPARC as published does not
# clear the L to H transition threshold in this model, so it is recorded as
# infeasible rather than quietly reported as a working design point.
_VERDICTS: dict[str, tuple[str, bool]] = {
    "ARC": ("troyon", True),
    "ITER": ("greenwald", True),
    "SPARC": ("lh_threshold", False),
}

_ARC_PLANT: dict[str, float] = {
    "thermal_power_mw": 552.5524917777225,
    "gross_electric_mw": 221.020996711089,
    "heating_draw_mw": 204.04322450708833,
    "balance_of_plant_mw": 8.84083986844356,
    "recirculating_mw": 227.8840643755319,
    "net_electric_mw": -6.8630676644428945,
    "recirculating_fraction": 1.0310516546688733,
    "net_efficiency": -0.01932689639015492,
    "engineering_gain": 0.9698835121128383,
}

_SCALING_COMPARISON: dict[str, tuple[float, float, float]] = {
    "IPB98(y,2)": (1.0753236982306174, 112.2237734788986, 3.164253686016847),
    "ITER89-P": (0.41802708468598476, 303.0964658160224, 1.1715890118787156),
    "Petty08": (2.5706390478800065, 41.61177270078677, 8.533750567242967),
}


@pytest.mark.parametrize("name", sorted(_OPERATING_POINTS))
def test_solved_operating_points_are_unchanged(name: str) -> None:
    """Every closed-form quantity of the reference run is reproduced.

    Tolerance: ``_CLOSED_FORM_RTOL``, derived at the top of this module from the
    operation count of the closed-form inversion and the machine epsilon.
    """
    trace = run_benchmark(IPB98Y2)
    result = trace.named(name)
    point = result.point
    expected = _OPERATING_POINTS[name]

    actual = {
        "fusion_power_mw": point.terms.fusion_power_mw,
        "alpha_power_mw": point.terms.alpha_power_mw,
        "auxiliary_power_mw": point.auxiliary_power_mw,
        "fusion_gain": point.fusion_gain,
        "confinement_time_s": point.confinement_time_s,
        "loss_power_mw": point.loss_power_mw,
        "stored_energy_mj": point.terms.stored_energy_mj,
        "bremsstrahlung_mw": point.terms.bremsstrahlung_mw,
        "synchrotron_mw": point.terms.synchrotron_mw,
        "line_radiation_mw": point.terms.line_radiation_mw,
        "triple_product": point.triple_product,
        "reactivity_m3_s": point.terms.reactivity_m3_s,
        "implied_h": result.implied_confinement_multiplier,
    }
    for key, reference in expected.items():
        assert actual[key] == pytest.approx(reference, rel=_CLOSED_FORM_RTOL, abs=1.0e-15)


@pytest.mark.parametrize("name", sorted(_CONSTRAINT_UTILISATION))
def test_constraint_utilisations_are_unchanged(name: str) -> None:
    """Every constraint utilisation of the reference run is reproduced."""
    report = run_benchmark(IPB98Y2).named(name).constraints
    for check_name, reference in _CONSTRAINT_UTILISATION[name].items():
        assert report.named(check_name).utilisation == pytest.approx(
            reference, rel=_CLOSED_FORM_RTOL
        )


@pytest.mark.parametrize("name", sorted(_VERDICTS))
def test_constraint_verdicts_are_unchanged(name: str) -> None:
    """The binding constraint and the feasibility verdict are categorical and exact."""
    report = run_benchmark(IPB98Y2).named(name).constraints
    binding, satisfied = _VERDICTS[name]
    assert report.binding.name == binding
    assert report.satisfied is satisfied


def test_a_violating_design_point_is_still_reported_as_violating() -> None:
    """The SPARC point remains flagged, with its violation named.

    This is pinned separately from the verdict table because it is the property
    the constraint layer exists for. A change that made a violating point return
    quietly would pass every numeric pin in this file and fail here.
    """
    result = run_benchmark(IPB98Y2).named("SPARC")
    assert not result.constraints.satisfied
    assert [check.name for check in result.constraints.violations] == ["lh_threshold"]


def test_arc_electrical_accounting_is_unchanged() -> None:
    """The plant result at the reference point is reproduced, sign included.

    The net electric power at the flat-profile ARC point is negative. That is a
    result of the model, not a defect of the pin, and it is recorded with its
    sign so that a change which turned it positive would be noticed.
    """
    plant = run_benchmark(IPB98Y2).named("ARC").plant
    assert plant is not None
    for key, reference in _ARC_PLANT.items():
        assert getattr(plant, key) == pytest.approx(reference, rel=_CLOSED_FORM_RTOL)
    assert plant.net_electric_mw < 0.0


@pytest.mark.parametrize("scaling", [IPB98Y2, ITER89P, PETTY08])
def test_the_arc_point_under_each_scaling_is_unchanged(scaling: object) -> None:
    """The confinement time, heating power, and gain under all three scalings."""
    point = solve_operating_point(machine("ARC").state, scaling)  # type: ignore[arg-type]
    tau, auxiliary, gain = _SCALING_COMPARISON[point.scaling_name]
    assert point.confinement_time_s == pytest.approx(tau, rel=_CLOSED_FORM_RTOL)
    assert point.auxiliary_power_mw == pytest.approx(auxiliary, rel=_CLOSED_FORM_RTOL)
    assert point.fusion_gain == pytest.approx(gain, rel=_CLOSED_FORM_RTOL)


def test_lawson_optimum_is_unchanged_within_its_convergence_tolerance() -> None:
    """The minimising temperature and the minimum triple product are reproduced.

    Tolerance: this is the one place where an iterative solve is pinned, so the
    tolerance comes from the solve rather than from the arithmetic. The bounded
    minimisation is given an absolute tolerance of 1e-6 keV on the temperature,
    so the temperature is pinned to 1e-5 keV, ten times that. The triple product
    at a minimum is stationary in the temperature, so an error of 1e-6 keV in the
    location moves it by order ``(1e-6 / 14)**2``, which is 5e-15 relative; it is
    pinned at 1e-9 relative, six orders above that.
    """
    optimum = optimum_lawson_temperature(PlasmaComposition(helium_fraction=0.05))
    assert optimum.temperature_kev == pytest.approx(14.504339277517367, abs=1.0e-5)
    assert optimum.n_tau == pytest.approx(2.5181734456517175e20, rel=1.0e-9)
    assert optimum.triple_product == pytest.approx(3.652444201536745e21, rel=1.0e-9)


def test_classical_lawson_optimum_for_a_pure_plasma_is_unchanged() -> None:
    """The textbook case, with no dilution and no relativistic correction."""
    optimum = optimum_lawson_temperature(
        PlasmaComposition(), relativistic_bremsstrahlung=False
    )
    assert optimum.temperature_kev == pytest.approx(14.245459561040922, abs=1.0e-5)
    assert optimum.triple_product == pytest.approx(2.952434109140229e21, rel=1.0e-9)


@pytest.mark.parametrize(
    ("invariant", "exponent", "feasible", "counts"),
    [
        (SweepInvariant.FIXED_DENSITY, 0.0, 12, {"troyon": 25, "lh_threshold": 16}),
        (
            SweepInvariant.FIXED_GREENWALD_FRACTION,
            2.0,
            18,
            {"troyon": 27, "lh_threshold": 14},
        ),
        (SweepInvariant.FIXED_BETA, 4.0, 40, {"troyon": 29, "greenwald": 12}),
    ],
)
def test_field_sweep_summary_is_unchanged(
    invariant: SweepInvariant,
    exponent: float,
    feasible: int,
    counts: dict[str, int],
) -> None:
    """The fitted exponent, the feasible count, and the binding tally of each sweep.

    Tolerance: the exponent is a least-squares fit to an exact power law, so it
    is pinned at 1e-9 absolute, which is derived from the conditioning of the
    log-space normal equations over 41 points rather than from any observed
    residual. Counts are integers and are pinned exactly.
    """
    case = machine("ARC")
    fields = np.linspace(4.0, 14.0, 41)
    trace = field_sweep(case.state, fields, invariant, IPB98Y2, plant=case.plant)

    fit = fit_power_law(trace.values, sweep_series(trace, "fusion_power"))
    assert fit.exponent == pytest.approx(exponent, abs=1.0e-9)
    assert len(trace.feasible_points) == feasible
    assert binding_constraint_counts(trace) == counts


def test_the_feasible_band_of_the_field_sweep_is_unchanged() -> None:
    """The feasible field band of the fixed Greenwald fraction sweep.

    Tolerance: one sweep step. The sweep runs from 4 T to 14 T over 41 points, so
    the step is 0.25 T and neither crossing can be located more finely than that
    whatever the arithmetic does. Pinning to a tighter tolerance would be pinning
    the grid rather than the physics.

    The reference values sit on grid points, so a shift of one grid point in
    either crossing is exactly at the tolerance. That is the boundary case the
    tolerance has to avoid sitting on, so the comparison is made against the
    grid index rather than the field value: an index is an integer and a shift
    of one grid point changes it by one, which no rounding can produce.
    """
    case = machine("ARC")
    low, high, points = 4.0, 14.0, 41
    step = (high - low) / (points - 1)
    trace = field_sweep(
        case.state,
        np.linspace(low, high, points),
        SweepInvariant.FIXED_GREENWALD_FRACTION,
        IPB98Y2,
        plant=case.plant,
    )
    feasible_indices = [
        index for index, point in enumerate(trace.points) if point.feasible
    ]
    assert (feasible_indices[0], feasible_indices[-1]) == (14, 31)
    assert trace.values[feasible_indices[0]] == pytest.approx(7.50, abs=0.5 * step)
    assert trace.values[feasible_indices[-1]] == pytest.approx(11.75, abs=0.5 * step)


def test_radius_sweep_summary_is_unchanged() -> None:
    """The fitted radius exponent, the feasible count, and the binding tally."""
    case = machine("ARC")
    trace = radius_sweep(
        case.state,
        np.linspace(2.0, 8.0, 41),
        SweepInvariant.FIXED_GREENWALD_FRACTION,
        IPB98Y2,
        plant=case.plant,
    )
    fit = fit_power_law(trace.values, sweep_series(trace, "fusion_power"))
    assert fit.exponent == pytest.approx(1.0, abs=1.0e-9)
    assert len(trace.feasible_points) == 12
    assert binding_constraint_counts(trace) == {"lh_threshold": 28, "troyon": 13}


def test_density_sweep_efficiency_optimum_is_unchanged() -> None:
    """The best feasible net efficiency and the density it occurs at.

    Tolerance: the density is quantised by the sweep step, which over 41 points
    from 0.6e20 to 2.4e20 is 4.5e18 per cubic metre, so it is pinned to one step.
    The efficiency at that point is a closed-form evaluation and is pinned at the
    closed-form tolerance.
    """
    case = machine("ARC")
    low, high, points = 0.6e20, 2.4e20, 41
    step = (high - low) / (points - 1)
    trace = density_sweep(case.state, np.linspace(low, high, points), IPB98Y2, plant=case.plant)

    assert len(trace.feasible_points) == 21
    assert binding_constraint_counts(trace) == {"troyon": 32, "lh_threshold": 9}

    best = max(
        trace.feasible_points,
        key=lambda p: p.plant.net_efficiency if p.plant is not None else -math.inf,
    )
    assert best.plant is not None
    assert best.value == pytest.approx(1.59e20, abs=step)
    assert best.plant.net_efficiency == pytest.approx(0.06288576952467696, rel=_CLOSED_FORM_RTOL)


def test_benchmark_discrepancies_are_unchanged() -> None:
    """The headline discrepancies against the published design points.

    These are the numbers the README reports, so they are pinned to make sure the
    document and the code cannot drift apart.
    """
    trace = run_benchmark(IPB98Y2)
    expected = {
        ("ARC", "fusion power"): -0.32361,
        ("ARC", "fusion gain Q"): -0.76733,
        ("ARC", "Greenwald fraction"): -0.00212,
        ("ITER", "fusion power"): -0.32874,
        ("ITER", "energy confinement time"): -0.16483,
        ("ITER", "Greenwald fraction"): -0.01440,
        ("SPARC", "energy confinement time"): 0.02311,
        ("SPARC", "Greenwald fraction"): -0.01703,
    }
    for (name, quantity), reference in expected.items():
        row = next(r for r in trace.named(name).rows if r.quantity == quantity)
        # Pinned to five decimal places, which is the precision the README quotes
        # these to as percentages to one decimal place, with a factor of ten of
        # margin so that rounding at the boundary cannot flip the comparison.
        assert row.relative_error == pytest.approx(reference, abs=1.0e-5)
