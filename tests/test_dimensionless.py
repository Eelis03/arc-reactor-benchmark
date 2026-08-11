"""Tier one: the dimensionless form of each confinement scaling.

Two kinds of check are applied. The first reaches outside this package: the ITER
Physics Basis quotes IPB98(y,2) in dimensionless variables as well as in
engineering ones, and reproducing the first set from the second tests the whole
conversion against numbers this repository did not produce.

The second is the perturbation check the rest of the suite uses on the scalings
themselves. Each dimensionless exponent is measured by scanning its own
parameter while the other three are held fixed, through the same operating point
solver every result in the package goes through, rather than by reading the
returned dataclass.
"""

from __future__ import annotations

import dataclasses
import itertools
import math

import pytest

from arc_benchmark.algorithm.balance import PlasmaState
from arc_benchmark.algorithm.dimensionless import dimensionless_exponents
from arc_benchmark.algorithm.operating import solve_operating_point
from arc_benchmark.model.confinement import (
    CONFINEMENT_SCALINGS,
    IPB98Y2,
    ITER89P,
    PETTY08,
    PowerLawConfinement,
)
from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.pipeline.machines import machine

_ARC = machine("ARC").state

# The dimensionless form of IPB98(y,2) as the ITER Physics Basis quotes it,
# alongside the engineering form this package transcribes. The two were produced
# by the same authors from the same regression, so recovering one from the other
# checks every exponent and the whole conversion at once.
_IPB98_PUBLISHED: dict[str, float] = {
    "normalised_gyroradius": -2.70,
    "beta": -0.90,
    "collisionality": -0.01,
    "safety_factor": -3.00,
}

_PUBLISHED_HALF_WIDTH = 0.005
"""Half-width of the rounding box around a published exponent.

Every engineering exponent in this package is quoted to two decimal places, so
the fit that produced it lies within this of the number printed.
"""

_ENGINEERING_NAMES = ("current", "field", "density", "power", "major_radius")

# Factors on the temperature, the field, the density, and the current that
# multiply one dimensionless parameter by f and leave the other three alone, at
# fixed geometry, shape, and isotope. Each row is the unique solution of the
# four definitions read as simultaneous conditions; the first is the standard
# rho* scan, in which a larger machine is imitated by a colder, weaker, and
# denser plasma.
_SCANS = {
    "normalised_gyroradius": lambda f: (f**-1.0, f**-1.5, f**-2.0, f**-1.5),
    "beta": lambda f: (f**0.5, f**0.25, f, f**0.25),
    "collisionality": lambda f: (f**-0.5, f**-0.25, 1.0, f**-0.25),
    "safety_factor": lambda f: (1.0, 1.0, 1.0, f**-1.0),
}

_SCAN_FACTOR = 1.2
"""Size of every scan below. Small enough to keep the temperature well inside
the range the reactivity fit is published over, and far enough from one that the
measured exponent is not a difference of nearly equal logarithms."""

# Tolerance on a measured exponent. The ratio below is a quotient of two
# closed-form solves of order fifty operations each, so it carries about 1e-15
# relative error, and dividing its logarithm by log(1.2) turns that into 6e-15
# on the exponent. The pin is 1e-12, two orders above that bound and derived
# from it rather than from any difference observed.
_MEASURED_ABS = 1.0e-12


def _rounding_spread(scaling: PowerLawConfinement, attribute: str) -> float:
    """How far one converted exponent moves over the rounding box of the published fit.

    Every tolerance in this file that compares against a published dimensionless
    exponent comes from here rather than from the difference that happens to be
    observed. The five engineering exponents are each quoted to two decimal
    places, so the underlying fit lies somewhere in a box of half-width 0.005
    about them, and walking the corners of that box bounds what the conversion
    can do with it.
    """
    reference = getattr(dimensionless_exponents(scaling), attribute)
    shifts = (-_PUBLISHED_HALF_WIDTH, 0.0, _PUBLISHED_HALF_WIDTH)
    worst = 0.0
    for corner in itertools.product(shifts, repeat=len(_ENGINEERING_NAMES)):
        moved = dataclasses.replace(
            scaling.exponents,
            **{
                name: getattr(scaling.exponents, name) + shift
                for name, shift in zip(_ENGINEERING_NAMES, corner, strict=True)
            },
        )
        converted = dimensionless_exponents(dataclasses.replace(scaling, exponents=moved))
        worst = max(worst, abs(getattr(converted, attribute) - reference))
    return worst


def _measured_exponent(scanned: PlasmaState, scaling: PowerLawConfinement) -> float:
    """Exponent of ``B tau_E`` between the reference ARC state and a scanned one."""
    reference = solve_operating_point(_ARC, scaling)
    solved = solve_operating_point(scanned, scaling)
    ratio = (scanned.toroidal_field * solved.confinement_time_s) / (
        _ARC.toroidal_field * reference.confinement_time_s
    )
    return math.log(ratio) / math.log(_SCAN_FACTOR)


@pytest.mark.parametrize("attribute", sorted(_IPB98_PUBLISHED))
def test_ipb98_reproduces_its_published_dimensionless_form(attribute: str) -> None:
    """The converted exponents match the ones the ITER Physics Basis prints.

    Tolerance: the width of the rounding box, computed by
    :func:`_rounding_spread`. The published dimensionless exponents came from the
    unrounded regression while this package carries the engineering exponents to
    two decimal places, so the two cannot agree more closely than that rounding
    allows. The bound is 0.11 on the gyroradius exponent and the observed
    difference is 0.007.
    """
    converted = getattr(dimensionless_exponents(IPB98Y2), attribute)
    assert converted == pytest.approx(
        _IPB98_PUBLISHED[attribute], abs=_rounding_spread(IPB98Y2, attribute)
    )


@pytest.mark.parametrize("scaling", sorted(CONFINEMENT_SCALINGS.values(), key=lambda s: s.name))
@pytest.mark.parametrize("attribute", sorted(_SCANS))
def test_each_exponent_is_recovered_by_scanning_its_own_parameter(
    scaling: PowerLawConfinement, attribute: str
) -> None:
    """Multiplying one dimensionless parameter by 1.2 moves ``B tau_E`` by its exponent.

    This is the perturbation check the confinement tests apply to the
    engineering exponents, applied to the dimensionless ones. It is stronger than
    reading the returned dataclass because it goes through the operating point
    solver, so it also asserts that the balance used to eliminate the loss power
    is the balance the solver actually inverts.

    Tolerance: ``_MEASURED_ABS``, derived at the top of this module.
    """
    temperature, field, density, current = _SCANS[attribute](_SCAN_FACTOR)
    scanned = dataclasses.replace(
        _ARC,
        temperature_kev=_ARC.temperature_kev * temperature,
        toroidal_field=_ARC.toroidal_field * field,
        electron_density=_ARC.electron_density * density,
        plasma_current_ma=_ARC.plasma_current_ma * current,
    )
    expected = getattr(dimensionless_exponents(scaling), attribute)
    assert _measured_exponent(scanned, scaling) == pytest.approx(expected, abs=_MEASURED_ABS)


@pytest.mark.parametrize("scaling", sorted(CONFINEMENT_SCALINGS.values(), key=lambda s: s.name))
def test_the_dimensional_residual_is_the_exponent_of_a_pure_size_scan(
    scaling: PowerLawConfinement,
) -> None:
    """Growing the machine at fixed dimensionless parameters moves ``B tau_E`` by the residual.

    This is what the Kadomtsev constraint says, stated as a measurement rather
    than as a definition. Scaling both radii by 1.2 while holding the normalised
    gyroradius, the beta, the collisionality, and the safety factor leaves a
    plasma that is dimensionlessly identical to the one it came from, so a
    scaling expressible in dimensionless variables alone would return the same
    Bohm-normalised confinement time. Whatever it returns instead is the
    residual, and the factors below are the unique solution of those four
    conditions at a size ratio of f.
    """
    geometry = _ARC.geometry
    scanned = dataclasses.replace(
        _ARC,
        geometry=PlasmaGeometry(
            major_radius=_SCAN_FACTOR * geometry.major_radius,
            minor_radius=_SCAN_FACTOR * geometry.minor_radius,
            elongation=geometry.elongation,
            triangularity=geometry.triangularity,
        ),
        temperature_kev=_ARC.temperature_kev * _SCAN_FACTOR**-0.5,
        toroidal_field=_ARC.toroidal_field * _SCAN_FACTOR**-1.25,
        electron_density=_ARC.electron_density * _SCAN_FACTOR**-2.0,
        plasma_current_ma=_ARC.plasma_current_ma * _SCAN_FACTOR**-0.25,
    )
    expected = dimensionless_exponents(scaling).dimensional_residual
    assert _measured_exponent(scanned, scaling) == pytest.approx(expected, abs=_MEASURED_ABS)


@pytest.mark.parametrize("scaling", [IPB98Y2, PETTY08])
def test_the_two_h_mode_fits_are_consistent_with_the_kadomtsev_constraint(
    scaling: PowerLawConfinement,
) -> None:
    """Their residuals are zero to within the rounding they are published at.

    IPB98(y,2) comes out at 0.008 against a bound of 0.086 and Petty08 at 0.019
    against 0.050, so both are consistent with satisfying the constraint exactly
    and this package cannot show otherwise from the two decimal places it
    carries. That is a statement about what the published form permits, not a
    claim that either fit satisfies it.
    """
    residual = dimensionless_exponents(scaling).dimensional_residual
    assert abs(residual) <= _rounding_spread(scaling, "dimensional_residual")


def test_iter89p_violates_the_kadomtsev_constraint_by_more_than_its_rounding() -> None:
    """The L-mode fit is the one that is genuinely not expressible in dimensionless form.

    Its residual is 0.075 against a rounding bound of 0.052, so unlike the other
    two it is not something the two decimal places of its published exponents can
    account for. The violation is small in size, and the second assertion is what
    it costs: doubling the machine at fixed normalised gyroradius, beta,
    collisionality, and safety factor moves the Bohm-normalised confinement time
    by 5.1 percent, where a scaling satisfying the constraint would leave it
    alone.
    """
    residual = dimensionless_exponents(ITER89P).dimensional_residual
    assert abs(residual) > _rounding_spread(ITER89P, "dimensional_residual")
    assert abs(2.0**residual - 1.0) == pytest.approx(0.051, abs=1.0e-3)


def test_petty_is_beta_independent_as_its_source_constrains_it_to_be() -> None:
    """Petty08 returns a beta exponent of zero, which is what it is fitted to be.

    This is the check the package could not previously make on its own
    description of that scaling. The measured exponent is 0.019 against a
    rounding bound of 0.031, so the residual is what two decimal places in the
    published exponents cost and not a failure of the constraint. The second
    assertion is the discriminating one: the beta dependence of IPB98(y,2) is
    almost fifty times larger, and it is that difference the constant-beta fit
    exists to isolate.
    """
    petty = dimensionless_exponents(PETTY08).beta
    assert abs(petty) <= _rounding_spread(PETTY08, "beta")
    assert abs(petty) < 0.1 * abs(dimensionless_exponents(IPB98Y2).beta)


def test_the_three_scalings_run_from_bohm_to_gyro_bohm() -> None:
    """ITER89-P is Bohm, Petty08 is gyro-Bohm, and IPB98(y,2) lies strictly between.

    Bohm diffusion gives a normalised gyroradius exponent of minus two and
    gyro-Bohm minus three, and the two outer scalings sit on those values to
    within the rounding of their published exponents. The bracketing of
    IPB98(y,2) needs no tolerance at all, which is why it is asserted as an
    inequality.
    """
    bohm = dimensionless_exponents(ITER89P).normalised_gyroradius
    gyro_bohm = dimensionless_exponents(PETTY08).normalised_gyroradius
    between = dimensionless_exponents(IPB98Y2).normalised_gyroradius

    assert abs(bohm + 2.0) <= _rounding_spread(ITER89P, "normalised_gyroradius")
    assert abs(gyro_bohm + 3.0) <= _rounding_spread(PETTY08, "normalised_gyroradius")
    assert gyro_bohm < between < bohm


def test_the_coefficient_does_not_enter_any_exponent() -> None:
    """Changing the coefficient moves the prefactor and nothing else."""
    rescaled = dataclasses.replace(IPB98Y2, coefficient=10.0 * IPB98Y2.coefficient)
    assert dimensionless_exponents(rescaled) == dimensionless_exponents(IPB98Y2)


@pytest.mark.parametrize("power", [0.0, -1.0, -1.5])
def test_a_scaling_with_no_steady_state_is_refused(power: float) -> None:
    """A power degradation outside ``(0, 1)`` leaves nothing to eliminate the loss power.

    The conversion divides by ``1 + a_P`` and the steady-state balance it comes
    from has no solution there, so the case is refused with the same message the
    operating point solver refuses it with rather than returning an exponent
    derived from a balance that does not close.
    """
    broken = dataclasses.replace(
        IPB98Y2, exponents=dataclasses.replace(IPB98Y2.exponents, power=power)
    )
    with pytest.raises(ValueError, match="power degradation"):
        dimensionless_exponents(broken)
