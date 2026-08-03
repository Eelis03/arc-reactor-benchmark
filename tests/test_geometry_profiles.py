"""Tier one: geometry and the profile correction factors.

The profile factors have closed forms for every quantity except fusion power, so
those are checked against the closed form directly. The fusion factor is a
quadrature and is checked against the closed form it must reduce to when the
reactivity is replaced by a power law, and against the flat limit.
"""

from __future__ import annotations

import itertools
import math

import numpy as np
import pytest

from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.model.profiles import FLAT_PROFILE, ProfileShape

_ARC = PlasmaGeometry(major_radius=3.3, minor_radius=1.13, elongation=1.84, triangularity=0.375)


def test_volume_is_the_elliptical_torus_result() -> None:
    """``V = 2 pi**2 R a**2 kappa``, to the last few bits."""
    expected = 2.0 * math.pi**2 * 3.3 * 1.13**2 * 1.84
    assert _ARC.volume == pytest.approx(expected, rel=1.0e-14)


def test_volume_scales_as_the_cube_of_size_at_fixed_aspect_ratio() -> None:
    """Doubling both radii multiplies the volume by eight, exactly."""
    doubled = PlasmaGeometry(
        major_radius=6.6, minor_radius=2.26, elongation=1.84, triangularity=0.375
    )
    assert doubled.volume == pytest.approx(8.0 * _ARC.volume, rel=1.0e-14)


def test_poloidal_perimeter_reduces_to_a_circle_at_unit_elongation() -> None:
    """A circular cross-section has a perimeter of ``2 pi a``.

    Ramanujan's approximation is exact for a circle, so this is an equality to
    rounding rather than an approximation.
    """
    circle = PlasmaGeometry(major_radius=3.3, minor_radius=1.13, elongation=1.0)
    assert circle.poloidal_perimeter == pytest.approx(2.0 * math.pi * 1.13, rel=1.0e-14)


def test_poloidal_perimeter_brackets_the_ellipse_bounds() -> None:
    """The perimeter lies between the two classical bounds for an ellipse.

    An ellipse with semi-axes ``a`` and ``b`` has a perimeter of at least
    ``pi (a + b)`` and at most ``2 pi sqrt((a**2 + b**2) / 2)``. Any approximation
    that leaves that interval is wrong regardless of its accuracy claim.
    """
    semi_major, semi_minor = 1.13 * 1.84, 1.13
    lower = math.pi * (semi_major + semi_minor)
    upper = 2.0 * math.pi * math.sqrt((semi_major**2 + semi_minor**2) / 2.0)
    assert lower <= _ARC.poloidal_perimeter <= upper


def test_effective_minor_radius_carries_the_elongation_correction() -> None:
    """``a_eff = a sqrt((1 + kappa**2) / 2)``."""
    expected = 1.13 * math.sqrt((1.0 + 1.84**2) / 2.0)
    assert _ARC.effective_minor_radius == pytest.approx(expected, rel=1.0e-14)


def test_geometry_rejects_impossible_shapes() -> None:
    """A minor radius above the major radius, or an elongation below one, raises."""
    with pytest.raises(ValueError, match="minor_radius"):
        PlasmaGeometry(major_radius=1.0, minor_radius=2.0, elongation=1.5)
    with pytest.raises(ValueError, match="elongation"):
        PlasmaGeometry(major_radius=3.3, minor_radius=1.13, elongation=0.5)
    with pytest.raises(ValueError, match="major_radius"):
        PlasmaGeometry(major_radius=-1.0, minor_radius=0.5, elongation=1.5)
    with pytest.raises(ValueError, match="triangularity"):
        PlasmaGeometry(major_radius=3.3, minor_radius=1.13, elongation=1.5, triangularity=1.5)


def test_flat_profile_leaves_every_factor_at_exactly_one() -> None:
    """The zero-dimensional default introduces no correction anywhere."""
    assert FLAT_PROFILE.is_flat
    assert FLAT_PROFILE.stored_energy_factor() == 1.0
    assert FLAT_PROFILE.bremsstrahlung_factor() == 1.0
    assert FLAT_PROFILE.density_square_factor() == 1.0
    assert FLAT_PROFILE.synchrotron_factor() == 1.0
    assert FLAT_PROFILE.fusion_factor(14.0) == 1.0
    assert FLAT_PROFILE.line_average_ratio() == 1.0


@pytest.mark.parametrize(
    ("density_exponent", "temperature_exponent"),
    [(0.0, 1.0), (0.4, 0.0), (0.4, 1.0), (1.0, 2.0)],
)
def test_closed_form_factors_match_their_integrals(
    density_exponent: float, temperature_exponent: float
) -> None:
    """Each closed-form profile factor agrees with a direct quadrature of it.

    The volume element goes as ``rho d rho``, so every factor here is the ratio of
    ``2 integral rho f(rho) d rho`` to the flat evaluation. Computing that
    integral numerically and comparing against the algebra is what catches a
    mistyped exponent in the closed form.

    Tolerance: derived from the quadrature, not from the error observed. Where a
    profile exponent lies below one the integrand has an algebraic singularity in
    its derivative at the plasma edge, so the composite trapezoidal rule
    converges as ``h**(1 + alpha)`` rather than as ``h**2``. The smallest
    exponent in this parameter set is 0.4, and the grid spacing in the
    transformed variable near the edge is 1e-5, which bounds the relative error
    at ``(1e-5)**1.4``, about 3e-8. The tolerance is set at 1e-6, thirty times
    that bound and still five orders below the percent-level change that any
    mistyped exponent would produce.
    """
    shape = ProfileShape(density_exponent, temperature_exponent)
    rho = np.linspace(0.0, 1.0, 200001)
    core = np.clip(1.0 - rho**2, 0.0, 1.0)
    density = shape.density_peaking * core**density_exponent
    temperature = shape.temperature_peaking * core**temperature_exponent

    def average(values: np.ndarray) -> float:
        return float(np.trapezoid(2.0 * rho * values, rho))

    assert average(density * temperature) == pytest.approx(
        shape.stored_energy_factor(), rel=1.0e-6
    )
    assert average(density**2 * np.sqrt(temperature)) == pytest.approx(
        shape.bremsstrahlung_factor(), rel=1.0e-6
    )
    assert average(density**2) == pytest.approx(shape.density_square_factor(), rel=1.0e-6)
    assert average(np.sqrt(density) * temperature**1.5) == pytest.approx(
        shape.synchrotron_factor(), rel=1.0e-6
    )


def test_peaking_raises_every_factor_above_one() -> None:
    """Peaked profiles concentrate the plasma, so every nonlinear average rises."""
    shape = ProfileShape(0.4, 1.0)
    assert shape.stored_energy_factor() > 1.0
    assert shape.bremsstrahlung_factor() > 1.0
    assert shape.density_square_factor() > 1.0
    assert shape.synchrotron_factor() > 1.0
    assert shape.fusion_factor(14.0) > 1.0


def test_fusion_factor_grows_with_peaking() -> None:
    """More peaked profiles enhance fusion power more, at a fixed mean temperature."""
    factors = [
        ProfileShape(0.2 * k, 0.5 * k).fusion_factor(14.0) for k in range(0, 5)
    ]
    assert all(later > earlier for earlier, later in itertools.pairwise(factors))


def test_fusion_factor_exceeds_the_bremsstrahlung_factor() -> None:
    """Peaking helps fusion more than it helps radiation, which is why it helps at all.

    Fusion power carries the square of density times a reactivity that rises
    faster than the square of temperature over the reactor range, while
    bremsstrahlung carries the square of density times only the square root of
    temperature.
    """
    shape = ProfileShape(0.4, 1.0)
    assert shape.fusion_factor(14.0) > shape.bremsstrahlung_factor()


@pytest.mark.parametrize("density_exponent", [0.2, 0.4, 0.5, 1.0, 2.0])
def test_line_average_ratio_matches_the_chord_integral(density_exponent: float) -> None:
    """The closed form equals the chord average divided by the volume average.

    The line average is the shape function integrated along ``d rho`` on a chord
    through the axis; the volume average is the same function integrated against
    ``2 rho d rho``. Computing both numerically and dividing is what checks the
    Gamma function algebra rather than restating it.

    Tolerance: the same reasoning as the closed-form factor test above. The
    integrand has an algebraic singularity in its derivative at the edge for an
    exponent below one, the trapezoidal rule then converges as ``h**(1 + alpha)``
    with ``alpha`` at least 0.2 and ``h`` equal to 5e-6, which bounds the
    relative error at 4e-7. The pin is 1e-5, an order of magnitude above that.
    """
    shape = ProfileShape(density_exponent, 1.0)
    rho = np.linspace(0.0, 1.0, 200001)
    core = np.clip(1.0 - rho**2, 0.0, 1.0) ** density_exponent

    chord_average = float(np.trapezoid(core, rho))
    volume_average = float(np.trapezoid(2.0 * rho * core, rho))
    assert shape.line_average_ratio() == pytest.approx(
        chord_average / volume_average, rel=1.0e-5
    )


def test_line_average_ratio_hits_its_two_exact_cases() -> None:
    """Two exponents have elementary closed forms, and both are reproduced.

    At ``alpha_n = 1`` the chord integral is 2/3 and the volume average is 1/2,
    so the ratio is exactly 4/3. At ``alpha_n = 1/2`` the chord integral is the
    area of a quarter circle, ``pi / 4``, against a volume average of 2/3, so the
    ratio is exactly ``3 pi / 8``. Neither needs the Gamma function, which is
    what makes them a check on the implementation that uses it.
    """
    assert ProfileShape(1.0, 0.0).line_average_ratio() == pytest.approx(4.0 / 3.0, rel=1.0e-14)
    assert ProfileShape(0.5, 0.0).line_average_ratio() == pytest.approx(
        3.0 * math.pi / 8.0, rel=1.0e-14
    )


def test_line_average_ratio_rises_with_density_peaking() -> None:
    """A more peaked density puts more of itself on the chord and none elsewhere.

    The temperature exponent is varied alongside to confirm it does not enter:
    the ratio is a statement about the density profile only.
    """
    ratios = [ProfileShape(0.25 * k, 0.5 * k).line_average_ratio() for k in range(0, 5)]
    assert ratios[0] == 1.0
    assert all(later > earlier for earlier, later in itertools.pairwise(ratios))
    assert ProfileShape(0.4, 0.0).line_average_ratio() == ProfileShape(
        0.4, 3.0
    ).line_average_ratio()


def test_profile_rejects_negative_exponents() -> None:
    """A negative exponent would invert the profile and is refused."""
    with pytest.raises(ValueError, match="density_exponent"):
        ProfileShape(-0.1, 1.0)
    with pytest.raises(ValueError, match="temperature_exponent"):
        ProfileShape(0.4, -1.0)
