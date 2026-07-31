"""Tier one: the operational limits, against published values for known machines.

The Greenwald limit, the Troyon limit, and the safety factor are all published
closed forms, so the strongest checks available are the ones that evaluate them
at a machine whose numbers are also published.
"""

from __future__ import annotations

import math

import pytest

from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.model.limits import (
    TROYON_COEFFICIENT,
    average_poloidal_field,
    bootstrap_fraction,
    current_at_fixed_q,
    cylindrical_safety_factor,
    greenwald_density,
    lh_threshold_power,
    normalised_beta,
    poloidal_beta,
    safety_factor_95,
    toroidal_beta,
    troyon_beta_limit,
)

_ITER_GEOMETRY = PlasmaGeometry(
    major_radius=6.2, minor_radius=2.0, elongation=1.70, triangularity=0.33
)
_ARC_GEOMETRY = PlasmaGeometry(
    major_radius=3.3, minor_radius=1.13, elongation=1.84, triangularity=0.375
)


def test_greenwald_limit_is_the_published_closed_form() -> None:
    """The limit equals ``I_p / (pi a**2)`` in units of 1e20 per cubic metre.

    Tolerance: this is one division and one multiplication, so agreement with the
    hand evaluation is to within a couple of units in the last place. 1e-14
    relative covers it.
    """
    computed = greenwald_density(15.0, 2.0)
    expected = 1.0e20 * 15.0 / (math.pi * 4.0)
    assert computed == pytest.approx(expected, rel=1.0e-14)
    assert computed == pytest.approx(1.1937e20, rel=1.0e-4)


@pytest.mark.parametrize(
    ("current_ma", "minor_radius", "density", "published_fraction"),
    [
        (15.0, 2.0, 1.0e20, 0.85),
        (7.8, 1.13, 1.3e20, 0.67),
        (8.7, 0.57, 3.1e20, 0.37),
    ],
)
def test_greenwald_fraction_matches_published_design_points(
    current_ma: float, minor_radius: float, density: float, published_fraction: float
) -> None:
    """ITER, ARC, and SPARC all report a Greenwald fraction this formula reproduces.

    Tolerance: the published fractions are quoted to two significant figures, so
    a value read from them carries a rounding uncertainty of half a unit in the
    second digit, which is between 0.7 and 1.4 percent depending on the value.
    The tolerance is set at 2 percent to cover the widest of those.
    """
    fraction = density / greenwald_density(current_ma, minor_radius)
    assert fraction == pytest.approx(published_fraction, rel=2.0e-2)


def test_troyon_limit_and_normalised_beta_are_exact_inverses() -> None:
    """Normalising the Troyon limit returns the limit coefficient itself.

    This round trip is what makes the two functions consistent: whatever the beta
    limit is, evaluating it and then normalising it must give back the normalised
    beta that defined it.
    """
    for beta_n_limit in (2.5, TROYON_COEFFICIENT, 3.5, 4.0):
        limit = troyon_beta_limit(2.0, 5.3, 15.0, beta_n_limit)
        recovered = normalised_beta(limit, 2.0, 5.3, 15.0)
        assert recovered == pytest.approx(beta_n_limit, rel=1.0e-14)


def test_toroidal_beta_is_the_pressure_ratio() -> None:
    """Beta equals twice the magnetic constant times pressure over field squared."""
    pressure, field = 5.0e5, 5.3
    expected = 2.0 * 1.25663706212e-6 * pressure / field**2
    assert toroidal_beta(pressure, field) == pytest.approx(expected, rel=1.0e-14)


def test_normalised_beta_of_iter_is_near_its_published_value() -> None:
    """A flat-profile ITER pressure gives a normalised beta near the published 1.8.

    Tolerance: 15 percent. The published value comes from a profile calculation
    with a fast alpha contribution that this volume-averaged pressure does not
    carry, so the two are not expected to agree closely. The assertion is that
    the formula lands in the right place, not that the model reproduces ITER.
    """
    density_sum = 1.0e20 + 0.9e20
    pressure = density_sum * 8.8 * 1.602176634e-16
    beta_t = toroidal_beta(pressure, 5.3)
    beta_n = normalised_beta(beta_t, 2.0, 5.3, 15.0)
    assert beta_n == pytest.approx(1.8, rel=0.15)


def test_safety_factor_95_reproduces_the_published_iter_value() -> None:
    """The shaping-corrected formula returns 3.0 at the ITER reference geometry.

    This is an external check on the shaping and aspect ratio correction factors,
    which are the parts of the formula most easily mistyped.

    Tolerance: the published ITER value is quoted as 3.0, to two significant
    figures, so the tolerance is half a unit in the second digit, 1.7 percent.
    The observed error is 0.13 percent.
    """
    q95 = safety_factor_95(_ITER_GEOMETRY, 5.3, 15.0)
    assert q95 == pytest.approx(3.0, rel=1.7e-2)


def test_cylindrical_safety_factor_and_its_inverse_round_trip() -> None:
    """Solving for the current at a chosen safety factor returns that factor."""
    for target_q in (2.0, 3.5, 5.0, 7.5):
        current = current_at_fixed_q(_ARC_GEOMETRY, 9.2, target_q)
        recovered = cylindrical_safety_factor(_ARC_GEOMETRY, 9.2, current)
        assert recovered == pytest.approx(target_q, rel=1.0e-13)


def test_safety_factor_scales_inversely_with_current() -> None:
    """Doubling the current halves the cylindrical safety factor, exactly."""
    base = cylindrical_safety_factor(_ARC_GEOMETRY, 9.2, 7.8)
    doubled = cylindrical_safety_factor(_ARC_GEOMETRY, 9.2, 15.6)
    assert doubled == pytest.approx(base / 2.0, rel=1.0e-14)


def test_poloidal_field_uses_the_true_perimeter_not_a_circle() -> None:
    """An elongated plasma has a longer boundary and therefore a weaker poloidal field.

    Using ``2 pi a`` in place of the true perimeter would overstate the poloidal
    field by the ratio of the two, which at an elongation of 1.84 is 45 percent,
    and would then overstate the bootstrap fraction by more than a factor of two.
    """
    circular_perimeter = 2.0 * math.pi * _ARC_GEOMETRY.minor_radius
    assert _ARC_GEOMETRY.poloidal_perimeter > circular_perimeter
    field = average_poloidal_field(_ARC_GEOMETRY, 7.8)
    naive = 1.25663706212e-6 * 7.8e6 / circular_perimeter
    assert field < naive
    assert field == pytest.approx(
        1.25663706212e-6 * 7.8e6 / _ARC_GEOMETRY.poloidal_perimeter, rel=1.0e-14
    )


def test_bootstrap_fraction_of_arc_is_near_its_published_value() -> None:
    """The ARC design point reports a bootstrap fraction near 0.63.

    Tolerance: 20 percent. The coefficient of the ``sqrt(eps) beta_p`` estimate
    is itself only known to within the range 0.5 to 1.0 for conventional
    profiles, which is a spread of plus or minus 36 percent about the 0.7 used
    here, so a tighter tolerance would be asserting more than the formula
    supports.
    """
    density_sum = 1.3e20 + 1.235e20
    pressure = density_sum * 14.0 * 1.602176634e-16
    beta_p = poloidal_beta(pressure, average_poloidal_field(_ARC_GEOMETRY, 7.8))
    assert bootstrap_fraction(_ARC_GEOMETRY, beta_p) == pytest.approx(0.63, rel=0.20)


def test_bootstrap_fraction_is_not_clipped_at_one() -> None:
    """A pressure gradient that would over-drive the current returns a value above one.

    Clipping would hide a design point that cannot exist, which is exactly the
    thing the constraint layer is there to report.
    """
    assert bootstrap_fraction(_ARC_GEOMETRY, beta_poloidal=5.0) > 1.0


def test_lh_threshold_matches_the_published_iter_estimate() -> None:
    """The Martin scaling gives an ITER threshold near the published 80 to 90 MW.

    Tolerance: the published ITER threshold is quoted across sources between 50
    and 90 MW depending on the density assumed at the transition, so the
    assertion is a range rather than a value. At the full baseline density of
    1e20 per cubic metre the scaling should land in the upper part of that range.
    """
    power = lh_threshold_power(1.0e20, 5.3, _ITER_GEOMETRY.surface_area)
    assert 70.0 < power < 100.0


def test_lh_threshold_scales_with_density_and_field_as_published() -> None:
    """The threshold carries the published exponents on density and field.

    Tolerance: two multiplications and two powers, so 1e-13 relative.
    """
    base = lh_threshold_power(1.0e20, 5.3, 683.0)
    denser = lh_threshold_power(2.0e20, 5.3, 683.0)
    stronger = lh_threshold_power(1.0e20, 10.6, 683.0)
    assert denser / base == pytest.approx(2.0**0.717, rel=1.0e-13)
    assert stronger / base == pytest.approx(2.0**0.803, rel=1.0e-13)


def test_limits_reject_unphysical_inputs() -> None:
    """Every limit refuses a non-positive current, field, or size."""
    with pytest.raises(ValueError, match="plasma_current_ma"):
        greenwald_density(0.0, 2.0)
    with pytest.raises(ValueError, match="minor_radius"):
        greenwald_density(15.0, 0.0)
    with pytest.raises(ValueError, match="toroidal_field"):
        toroidal_beta(1.0e5, 0.0)
    with pytest.raises(ValueError, match="beta_n_limit"):
        troyon_beta_limit(2.0, 5.3, 15.0, 0.0)
    with pytest.raises(ValueError, match="cylindrical_q"):
        current_at_fixed_q(_ARC_GEOMETRY, 9.2, 0.0)
