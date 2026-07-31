"""Tier one: the Lawson condition, ignition, and the triple product."""

from __future__ import annotations

import itertools
import math

import pytest

from arc_benchmark.algorithm.balance import PlasmaComposition
from arc_benchmark.algorithm.lawson import (
    lawson_n_tau,
    lawson_triple_product,
    optimum_lawson_temperature,
)
from arc_benchmark.model.constants import ALPHA_FRACTION
from arc_benchmark.model.radiation import ImpurityRadiator

_PURE = PlasmaComposition()
_ASHY = PlasmaComposition(helium_fraction=0.05)


def test_minimum_triple_product_matches_the_classical_value() -> None:
    """A pure deuterium tritium plasma needs about 3e21 m^-3 keV s to ignite.

    The classical result, quoted in every fusion textbook, is a minimum ignition
    triple product of roughly 3e21 inverse cubic metre keV second at a
    temperature between 15 and 30 keV. Deriving it here from the same power
    balance the rest of the package solves, rather than quoting it, is what makes
    it a check.

    Tolerance: the textbook value is quoted to one significant figure and the
    temperature to a range, so the triple product is checked to 25 percent and
    the temperature to the interval 10 to 30 keV. A tighter tolerance would be
    asserting more precision than the reference carries.
    """
    optimum = optimum_lawson_temperature(_PURE, relativistic_bremsstrahlung=False)
    assert optimum.triple_product == pytest.approx(3.0e21, rel=0.25)
    assert 10.0 < optimum.temperature_kev < 30.0


def test_ignition_requires_more_than_a_gain_of_ten() -> None:
    """The requirement rises monotonically with the target gain, up to ignition."""
    values = [
        lawson_n_tau(15.0, _PURE, gain) for gain in (1.0, 2.0, 5.0, 10.0, 50.0, math.inf)
    ]
    assert all(b > a for a, b in itertools.pairwise(values))


def test_the_gain_multiplier_enters_exactly_as_derived() -> None:
    """At a finite gain the heating term carries ``1 + 1 / (f_alpha Q)``.

    The requirement is ``S / (g_Q A - R)`` with ``S`` the stored energy
    coefficient, ``A`` the alpha heating coefficient, ``R`` the radiation
    coefficient, and ``g_Q = 1 + 1 / (f_alpha Q)``. Taking the reciprocal of the
    requirement at a finite gain and subtracting the reciprocal at ignition
    cancels ``R`` entirely and leaves ``A / (S f_alpha Q)``, so the same
    difference taken at two gains must be in the inverse ratio of those gains.
    Nothing in that identity depends on any coefficient, which is what makes it a
    check on the algebra rather than a restatement of it.

    Tolerance: the reconstruction is four reciprocals, two subtractions, and one
    division on numbers of order 1e-20, so 1e-11 relative covers the
    cancellation with several orders to spare.
    """
    temperature = 15.0
    at_ignition = 1.0 / lawson_n_tau(temperature, _PURE, math.inf)
    at_two = 1.0 / lawson_n_tau(temperature, _PURE, 2.0)
    at_eight = 1.0 / lawson_n_tau(temperature, _PURE, 8.0)

    ratio = (at_two - at_ignition) / (at_eight - at_ignition)
    assert ratio == pytest.approx(8.0 / 2.0, rel=1.0e-11)

    # The alpha fraction is 0.19999, not exactly a fifth, because the reaction
    # releases 3.518 MeV to the alpha out of 17.589 MeV in total.
    assert pytest.approx(0.2, rel=1.0e-4) == ALPHA_FRACTION


def test_dilution_raises_the_requirement() -> None:
    """Helium ash both removes fuel and raises the effective charge, so it costs twice."""
    assert lawson_n_tau(15.0, _ASHY) > lawson_n_tau(15.0, _PURE)


def test_a_cold_plasma_cannot_ignite_at_any_confinement() -> None:
    """Below the ideal ignition temperature the requirement is infinite.

    At low temperature the bremsstrahlung loss exceeds the alpha heating whatever
    the density, so no confinement time is sufficient. Returning infinity is the
    correct answer; returning a large finite number would not be.
    """
    assert math.isinf(lawson_n_tau(2.0, _PURE))
    assert math.isfinite(lawson_n_tau(10.0, _PURE))


def test_an_impurity_load_can_make_ignition_impossible() -> None:
    """Enough line radiation removes the ignition window entirely."""
    poisoned = PlasmaComposition(
        impurities=(ImpurityRadiator("tungsten", 74, 5.0e-4, 8.0e-32),)
    )
    assert math.isinf(lawson_n_tau(15.0, poisoned))
    assert math.isfinite(lawson_n_tau(15.0, PlasmaComposition()))


def test_line_radiation_can_be_excluded_from_the_condition() -> None:
    """Turning the line term off lowers the requirement, since it is a loss."""
    dirty = PlasmaComposition(
        impurities=(ImpurityRadiator("argon", 18, 1.0e-3, 3.0e-33),)
    )
    with_line = lawson_n_tau(15.0, dirty, include_line_radiation=True)
    without_line = lawson_n_tau(15.0, dirty, include_line_radiation=False)
    assert without_line < with_line


def test_triple_product_is_the_requirement_times_the_temperature() -> None:
    """The two forms of the condition are consistent, to the last bit."""
    for temperature in (8.0, 15.0, 25.0, 40.0):
        point = lawson_triple_product(temperature, _ASHY)
        assert point.triple_product == pytest.approx(
            point.n_tau * temperature, rel=1.0e-15
        )


def test_the_optimum_is_a_genuine_minimum() -> None:
    """The returned temperature beats both of its neighbours.

    Tolerance: the minimiser is given an absolute tolerance of 1e-6 keV on the
    temperature, so the neighbours are taken one keV away, six orders outside
    that tolerance, and the comparison is a strict inequality.
    """
    optimum = optimum_lawson_temperature(_ASHY)
    for offset in (-1.0, 1.0):
        neighbour = lawson_triple_product(optimum.temperature_kev + offset, _ASHY)
        assert neighbour.triple_product > optimum.triple_product


def test_relativistic_correction_raises_the_requirement() -> None:
    """The correction adds a loss, so the plasma has to work harder to ignite."""
    classical = lawson_n_tau(25.0, _PURE, relativistic_bremsstrahlung=False)
    corrected = lawson_n_tau(25.0, _PURE, relativistic_bremsstrahlung=True)
    assert corrected > classical


def test_lawson_rejects_a_non_positive_gain() -> None:
    """A gain of zero or below is not a target the condition can be written for."""
    with pytest.raises(ValueError, match="gain"):
        lawson_n_tau(15.0, _PURE, 0.0)


def test_optimum_rejects_a_bracket_outside_the_reactivity_range() -> None:
    """A bracket that does not intersect the published fit range raises."""
    with pytest.raises(ValueError, match="bracket"):
        optimum_lawson_temperature(_PURE, bracket_kev=(200.0, 400.0))
