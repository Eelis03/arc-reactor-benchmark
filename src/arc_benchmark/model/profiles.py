"""Profile correction factors for a volume-averaged model.

A zero-dimensional balance evaluates every term at the volume-averaged density
and temperature. The true volume integrals are not equal to those evaluations
because fusion power goes as the square of density times a steeply rising
function of temperature, and both quantities are peaked on axis. The ratio
between the integral and the flat-profile evaluation is a pure number that
depends only on the profile shape, so it can be computed once and applied as a
multiplier without giving up the zero-dimensional structure.

Profiles are the standard parabolic family

    n(rho) = n0 (1 - rho**2)**alpha_n,  T(rho) = T0 (1 - rho**2)**alpha_T

with ``rho`` the normalised minor radius. For the elliptical cross-section used
throughout, the volume element is proportional to ``rho d rho``, so the volume
average of ``(1 - rho**2)**alpha`` is ``1 / (1 + alpha)``.

One factor here is not a volume integral at all.
:meth:`ProfileShape.line_average_ratio` converts the volume-averaged density
this model carries into the line-averaged density that the confinement scalings,
the Greenwald limit, and the L to H threshold are all published against. It is a
ratio of two integrals of the same shape function over different measures, so it
is the same kind of pure number and is applied in the same way.

Setting both exponents to zero recovers flat profiles and every factor here
becomes exactly one, which is how the rest of the package stays honest about
what a genuinely zero-dimensional answer looks like.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np
from scipy import integrate

from arc_benchmark.model.reactivity import REACTIVITY_TEMPERATURE_RANGE_KEV, dt_reactivity

__all__ = ["FLAT_PROFILE", "ProfileShape"]

_QUADRATURE_POINTS: Final[int] = 801
"""Odd, so the composite Simpson rule applies. Fixed rather than adaptive so
that the same grid, and therefore the same answer, is used on every machine."""


@dataclass(frozen=True, slots=True)
class ProfileShape:
    """Parabolic density and temperature profile exponents.

    Attributes:
        density_exponent: ``alpha_n``. Zero gives a flat density.
        temperature_exponent: ``alpha_T``. Zero gives a flat temperature.
    """

    density_exponent: float = 0.0
    temperature_exponent: float = 0.0

    def __post_init__(self) -> None:
        """Reject exponents that would invert the profile."""
        if self.density_exponent < 0.0:
            raise ValueError(f"density_exponent must be non-negative, got {self.density_exponent}")
        if self.temperature_exponent < 0.0:
            raise ValueError(
                f"temperature_exponent must be non-negative, got {self.temperature_exponent}"
            )

    @property
    def is_flat(self) -> bool:
        """True when both exponents are zero and every factor here is one."""
        return self.density_exponent == 0.0 and self.temperature_exponent == 0.0

    @property
    def density_peaking(self) -> float:
        """``n0 / <n>``, equal to ``1 + alpha_n``."""
        return 1.0 + self.density_exponent

    @property
    def temperature_peaking(self) -> float:
        """``T0 / <T>``, equal to ``1 + alpha_T``."""
        return 1.0 + self.temperature_exponent

    def stored_energy_factor(self) -> float:
        """``<n T> / (<n> <T>)``, in closed form.

        The integral of ``(1 - rho**2)**(alpha_n + alpha_T)`` gives
        ``(1 + alpha_n)(1 + alpha_T) / (1 + alpha_n + alpha_T)``. Peaked profiles
        store more energy than the product of the averages suggests, so this is
        always at least one.
        """
        return (
            self.density_peaking
            * self.temperature_peaking
            / (1.0 + self.density_exponent + self.temperature_exponent)
        )

    def density_square_factor(self) -> float:
        """``<n**2> / <n>**2``, in closed form.

        Equal to ``(1 + alpha_n)**2 / (1 + 2 alpha_n)``. This is the factor that
        applies to any loss going as the square of density with no temperature
        dependence, which in this package means impurity line radiation.
        """
        return self.density_peaking**2 / (1.0 + 2.0 * self.density_exponent)

    def synchrotron_factor(self) -> float:
        """``<n**0.5 T**1.5> / (<n>**0.5 <T>**1.5)``, in closed form.

        The synchrotron loss after reabsorption carries a half power of density
        and a one and a half power of temperature, so peaking raises it. The
        field enters only as a constant multiplier and drops out of the ratio.
        """
        exponent = 0.5 * self.density_exponent + 1.5 * self.temperature_exponent
        return float(
            math.sqrt(self.density_peaking) * self.temperature_peaking**1.5 / (1.0 + exponent)
        )

    def bremsstrahlung_factor(self) -> float:
        """``<n**2 sqrt(T)> / (<n>**2 sqrt(<T>))``, in closed form.

        The exponent of ``(1 - rho**2)`` in the integrand is
        ``2 alpha_n + alpha_T / 2``.
        """
        exponent = 2.0 * self.density_exponent + 0.5 * self.temperature_exponent
        return self.density_peaking**2 * math.sqrt(self.temperature_peaking) / (1.0 + exponent)

    def line_average_ratio(self) -> float:
        """``n_line / <n>``, the line average over the volume average, in closed form.

        Every confinement scaling in this package, the Greenwald limit, and the
        L to H power threshold are published against the density an
        interferometer measures along a chord, not against the volume average a
        zero-dimensional model carries. For the parabolic family both averages
        are integrals of the same shape function, along ``d rho`` on the chord
        and against ``2 rho d rho`` over the volume, so their ratio is a pure
        number that depends only on ``alpha_n``:

            n_line / <n> = (1 + alpha_n) sqrt(pi) Gamma(1 + alpha_n)
                           / (2 Gamma(alpha_n + 3/2))

        from the Beta function identity
        ``int_0^1 (1 - x**2)**alpha dx = sqrt(pi) Gamma(1 + alpha) / (2 Gamma(alpha + 3/2))``.
        It is 1.0 for a flat density, 1.1447 at ``alpha_n = 0.4``, and 1.1781 at
        ``alpha_n = 0.5``.

        The chord is taken to pass through the magnetic axis of nested elliptical
        flux surfaces of constant elongation, which makes the chord coordinate
        proportional to ``rho`` whether the chord is vertical or on the midplane.
        A Shafranov-shifted equilibrium, or a tangential chord, would break that
        proportionality and is not modelled.

        Returns:
            The ratio, at least one, and exactly one for a flat density.
        """
        if self.density_exponent == 0.0:
            return 1.0
        return float(
            self.density_peaking
            * math.sqrt(math.pi)
            * math.gamma(1.0 + self.density_exponent)
            / (2.0 * math.gamma(self.density_exponent + 1.5))
        )

    def fusion_factor(self, mean_temperature_kev: float) -> float:
        """``<n**2 sigma_v(T)> / (<n>**2 sigma_v(<T>))``, by quadrature.

        The reactivity is not a power law, so this factor depends on the mean
        temperature as well as on the profile shape, and no closed form exists.
        A composite Simpson rule on a fixed grid of
        :data:`_QUADRATURE_POINTS` nodes is used.

        Where the local temperature falls below the lower end of the published
        Bosch and Hale range the reactivity is set to zero rather than
        extrapolated. At 0.2 keV the reactivity is nine orders of magnitude below
        its value at 14 keV and the density weighting there is also small, so the
        truncation changes the factor by far less than the quadrature error.

        Args:
            mean_temperature_kev: Volume-averaged temperature in keV.

        Returns:
            The dimensionless enhancement, exactly one for a flat profile.
        """
        if self.is_flat:
            return 1.0

        rho = np.linspace(0.0, 1.0, _QUADRATURE_POINTS, dtype=np.float64)
        shape = np.clip(1.0 - rho**2, 0.0, 1.0)

        density = shape**self.density_exponent
        temperature = (
            self.temperature_peaking * mean_temperature_kev * shape ** (self.temperature_exponent)
        )

        low, high = REACTIVITY_TEMPERATURE_RANGE_KEV
        inside = (temperature >= low) & (temperature <= high)
        sigma_v = np.zeros_like(temperature)
        if np.any(inside):
            sigma_v[inside] = dt_reactivity(temperature[inside])

        weighted = 2.0 * rho * (self.density_peaking * density) ** 2 * sigma_v
        integral = float(integrate.simpson(weighted, x=rho))
        reference = float(dt_reactivity(mean_temperature_kev))
        return integral / reference


FLAT_PROFILE: Final[ProfileShape] = ProfileShape()
"""The zero-dimensional default: every correction factor is exactly one."""
