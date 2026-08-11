"""Tier one: the radiation terms, their coefficients, and their scalings."""

from __future__ import annotations

import math

import pytest

from arc_benchmark.model.radiation import (
    BREMSSTRAHLUNG_COEFFICIENT,
    LARMOR_COEFFICIENT,
    ImpurityRadiator,
    bremsstrahlung_density,
    cyclotron_escape_factor,
    line_radiation_density,
    optically_thin_cyclotron_density,
    synchrotron_density,
)


def test_classical_bremsstrahlung_is_the_hand_evaluation() -> None:
    """The classical form reproduces an arithmetic evaluation exactly.

    Tolerance: the computation is three multiplications and one square root, so
    the result should agree with a hand evaluation to within a few units in the
    last place. 1e-14 relative covers that with two orders of magnitude to spare.
    """
    n_e, t_e, z_eff = 1.0e20, 9.0, 1.5
    expected = BREMSSTRAHLUNG_COEFFICIENT * z_eff * n_e**2 * math.sqrt(t_e)
    computed = float(bremsstrahlung_density(n_e, t_e, z_eff, relativistic=False))
    assert computed == pytest.approx(expected, rel=1.0e-14)


def test_bremsstrahlung_scales_as_the_square_of_density() -> None:
    """Doubling density quadruples bremsstrahlung, exactly."""
    base = float(bremsstrahlung_density(1.0e20, 10.0, 1.0, relativistic=False))
    doubled = float(bremsstrahlung_density(2.0e20, 10.0, 1.0, relativistic=False))
    assert doubled == pytest.approx(4.0 * base, rel=1.0e-14)


def test_bremsstrahlung_scales_as_the_square_root_of_temperature() -> None:
    """Quadrupling temperature doubles bremsstrahlung, exactly."""
    base = float(bremsstrahlung_density(1.0e20, 5.0, 1.0, relativistic=False))
    hotter = float(bremsstrahlung_density(1.0e20, 20.0, 1.0, relativistic=False))
    assert hotter == pytest.approx(2.0 * base, rel=1.0e-14)


def test_relativistic_correction_is_small_and_positive() -> None:
    """The correction raises the loss by about 6 percent at 14 keV.

    Tolerance: the multiplier is evaluated in closed form, so the assertion is on
    its size rather than on a fitted number. It is checked to lie between 5 and
    8 percent, which is wide enough to be a statement about the physics and
    narrow enough to catch a wrong power of the electron rest energy.
    """
    classical = float(bremsstrahlung_density(1.0e20, 14.0, 1.0, relativistic=False))
    corrected = float(bremsstrahlung_density(1.0e20, 14.0, 1.0, relativistic=True))
    ratio = corrected / classical
    assert 1.05 < ratio < 1.08


def test_bremsstrahlung_rejects_an_effective_charge_below_one() -> None:
    """An effective charge below one is unphysical and raises."""
    with pytest.raises(ValueError, match="z_effective"):
        bremsstrahlung_density(1.0e20, 10.0, 0.5)


def test_larmor_coefficient_matches_the_literature_value() -> None:
    """The assembled constant reproduces the 6.2e-3 usually quoted.

    The optically thin cyclotron loss is often written as
    ``6.2e-3 B**2 n_20 T_keV`` in megawatt per cubic metre. Assembling the same
    coefficient from the elementary charge, the electron mass, the electric
    constant, and the speed of light must give the same number.

    Tolerance: the quoted literal carries two significant figures, so it is
    checked to 1 percent.
    """
    quoted_in_si = 6.2e-3 * 1.0e6 / 1.0e20
    assert pytest.approx(quoted_in_si, rel=1.0e-2) == LARMOR_COEFFICIENT


def test_optically_thin_cyclotron_scales_as_the_square_of_field() -> None:
    """Doubling field quadruples the unreabsorbed cyclotron emission, exactly."""
    base = float(optically_thin_cyclotron_density(1.0e20, 10.0, 5.0))
    doubled = float(optically_thin_cyclotron_density(1.0e20, 10.0, 10.0))
    assert doubled == pytest.approx(4.0 * base, rel=1.0e-14)


def test_escape_factor_is_a_small_fraction_at_reactor_parameters() -> None:
    """A reactor core reabsorbs the overwhelming majority of its cyclotron emission.

    Tolerance: the assertion is a range rather than a value, because the point is
    that the factor is of order a few percent. A factor at or above one would
    mean the formula had been written as a loss rather than as a fraction.
    """
    for density, temperature, field, minor_radius in (
        (1.0e20, 8.8, 5.3, 2.00),
        (1.3e20, 14.0, 9.2, 1.13),
        (3.1e20, 7.3, 12.2, 0.57),
    ):
        factor = float(cyclotron_escape_factor(density, temperature, field, minor_radius, 0.9))
        assert 0.005 < factor < 0.1


def test_synchrotron_is_the_product_of_the_two_parts() -> None:
    """The exported loss is exactly the thin emission times the escaping fraction."""
    args = (1.3e20, 14.0, 9.2)
    thin = float(optically_thin_cyclotron_density(*args))
    factor = float(cyclotron_escape_factor(*args, 1.13, 0.9))
    assert float(synchrotron_density(*args, 1.13, 0.9)) == pytest.approx(thin * factor, rel=1.0e-14)


def test_synchrotron_carries_the_expected_field_exponent() -> None:
    """Doubling the field multiplies the synchrotron loss by two to the five halves.

    The thin emission carries two powers of field and the escape factor carries a
    further half power, so the product carries two and a half. This is the
    exponent that makes synchrotron the cost of the high-field pathway, so it is
    checked directly rather than inferred.
    """
    low = float(synchrotron_density(1.3e20, 14.0, 5.0, 1.13, 0.9))
    high = float(synchrotron_density(1.3e20, 14.0, 10.0, 1.13, 0.9))
    assert high / low == pytest.approx(2.0**2.5, rel=1.0e-13)


def test_synchrotron_carries_the_expected_temperature_exponent() -> None:
    """Doubling the temperature multiplies the loss by two to the three halves."""
    low = float(synchrotron_density(1.3e20, 7.0, 9.2, 1.13, 0.9))
    high = float(synchrotron_density(1.3e20, 14.0, 9.2, 1.13, 0.9))
    assert high / low == pytest.approx(2.0**1.5, rel=1.0e-13)


def test_synchrotron_carries_the_expected_density_exponent() -> None:
    """Doubling the density multiplies the loss by the square root of two."""
    low = float(synchrotron_density(1.0e20, 14.0, 9.2, 1.13, 0.9))
    high = float(synchrotron_density(2.0e20, 14.0, 9.2, 1.13, 0.9))
    assert high / low == pytest.approx(math.sqrt(2.0), rel=1.0e-13)


def test_reflectivity_enters_as_the_square_root_of_what_is_not_reflected() -> None:
    """A wall reflecting three quarters halves the loss, and one reflecting all removes it.

    The escape factor carries ``sqrt(1 - R_w)``, so the ratio between two
    reflectivities is the square root of the ratio of what each lets through.

    Both reflectivities are chosen so that ``1 - R_w`` is exact in binary
    floating point: 0.75 is exact outright, and the near-perfect wall is placed
    at ``1 - 2**-20`` rather than at a decimal value near one. A decimal such as
    ``1 - 1e-12`` is not representable, and the subtraction then loses four
    significant digits to cancellation, so a test written that way would be
    measuring the cancellation rather than the square root.

    Tolerance: with the subtraction exact, the remaining arithmetic is one square
    root and one division, so 1e-13 relative.
    """
    none_reflected = float(synchrotron_density(1.3e20, 14.0, 9.2, 1.13, 0.0))
    three_quarters = float(synchrotron_density(1.3e20, 14.0, 9.2, 1.13, 0.75))
    almost_all = float(synchrotron_density(1.3e20, 14.0, 9.2, 1.13, 1.0 - 2.0**-20))

    assert none_reflected > 0.0
    assert three_quarters / none_reflected == pytest.approx(0.5, rel=1.0e-13)
    assert almost_all / none_reflected == pytest.approx(2.0**-10, rel=1.0e-13)


def test_synchrotron_rejects_an_impossible_reflectivity() -> None:
    """A reflectivity at or above one, or below zero, raises."""
    with pytest.raises(ValueError, match="wall_reflectivity"):
        synchrotron_density(1.0e20, 10.0, 5.0, 2.0, 1.0)
    with pytest.raises(ValueError, match="wall_reflectivity"):
        synchrotron_density(1.0e20, 10.0, 5.0, 2.0, -0.1)


def test_line_radiation_sums_over_species() -> None:
    """Two species radiate the sum of what each would radiate alone."""
    argon = ImpurityRadiator("argon", 18, 0.001, 3.0e-33)
    tungsten = ImpurityRadiator("tungsten", 74, 1.0e-5, 8.0e-32)
    n_e = 1.0e20
    separately = float(line_radiation_density(n_e, (argon,))) + float(
        line_radiation_density(n_e, (tungsten,))
    )
    together = float(line_radiation_density(n_e, (argon, tungsten)))
    assert together == pytest.approx(separately, rel=1.0e-14)


def test_line_radiation_of_a_clean_plasma_is_zero() -> None:
    """No impurities means no line radiation."""
    assert float(line_radiation_density(1.0e20, ())) == 0.0


def test_impurity_rejects_unphysical_parameters() -> None:
    """Charge, concentration, and loss parameter are all validated."""
    with pytest.raises(ValueError, match="atomic_number"):
        ImpurityRadiator("nothing", 0, 0.01, 1.0e-32)
    with pytest.raises(ValueError, match="concentration"):
        ImpurityRadiator("argon", 18, 1.5, 1.0e-32)
    with pytest.raises(ValueError, match="loss_parameter"):
        ImpurityRadiator("argon", 18, 0.01, -1.0)
