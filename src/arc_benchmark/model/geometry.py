"""Plasma geometry: the shape quantities every other model term is built on.

A zero-dimensional model needs three geometric numbers, the plasma volume, the
poloidal circumference, and the inverse aspect ratio, plus the elongation that
the confinement scalings and the safety factor both depend on. They are collected
here so that no other module computes a volume of its own.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

__all__ = ["PlasmaGeometry"]


@dataclass(frozen=True, slots=True)
class PlasmaGeometry:
    """An elongated, D-shaped plasma cross-section.

    Attributes:
        major_radius: Geometric major radius of the magnetic axis in metre.
        minor_radius: Horizontal minor radius in metre.
        elongation: Ratio of vertical to horizontal half-height, dimensionless.
        triangularity: Plasma triangularity, dimensionless.
    """

    major_radius: float
    minor_radius: float
    elongation: float
    triangularity: float = 0.0

    def __post_init__(self) -> None:
        """Reject geometries that are not physical."""
        if self.major_radius <= 0.0:
            raise ValueError(f"major_radius must be positive, got {self.major_radius}")
        if not 0.0 < self.minor_radius < self.major_radius:
            raise ValueError(
                "minor_radius must be positive and smaller than major_radius, got "
                f"{self.minor_radius} against {self.major_radius}"
            )
        if self.elongation < 1.0:
            raise ValueError(f"elongation must be at least 1, got {self.elongation}")
        if not -1.0 < self.triangularity < 1.0:
            raise ValueError(f"triangularity must lie in (-1, 1), got {self.triangularity}")

    @property
    def inverse_aspect_ratio(self) -> float:
        """``a / R``, written ``epsilon`` in the confinement scalings."""
        return self.minor_radius / self.major_radius

    @property
    def aspect_ratio(self) -> float:
        """``R / a``."""
        return self.major_radius / self.minor_radius

    @property
    def volume(self) -> float:
        """Plasma volume in cubic metre.

        The elliptical torus result ``2 pi**2 R a**2 kappa``. Triangularity moves
        this by well under one percent for the shapes considered here and is not
        included, which keeps the volume an exact function of three numbers.
        """
        return 2.0 * math.pi**2 * self.major_radius * self.minor_radius**2 * self.elongation

    @property
    def surface_area(self) -> float:
        """Plasma surface area in square metre.

        Uses the Ramanujan approximation for the poloidal perimeter of the
        elliptical cross-section, swept about the major axis.
        """
        return 2.0 * math.pi * self.major_radius * self.poloidal_perimeter

    @property
    def poloidal_perimeter(self) -> float:
        """Poloidal circumference in metre, Ramanujan's second approximation.

        For an ellipse with semi-axes ``a`` and ``kappa a``. The approximation is
        accurate to better than 1e-5 relative for every elongation below 3.
        """
        semi_major = self.minor_radius * self.elongation
        semi_minor = self.minor_radius
        h = ((semi_major - semi_minor) / (semi_major + semi_minor)) ** 2
        return (
            math.pi
            * (semi_major + semi_minor)
            * (1.0 + 3.0 * h / (10.0 + math.sqrt(4.0 - 3.0 * h)))
        )

    @property
    def areal_elongation(self) -> float:
        """``kappa_a = S / (pi a**2)``, the elongation the IPB98 fit uses.

        For the elliptical cross-section assumed here the poloidal cross-section
        area is ``pi a**2 kappa``, so this equals :attr:`elongation`. It is kept
        as a separate name because the published scaling is written in terms of
        the areal quantity and the two differ for a shaped boundary.
        """
        return self.elongation

    @property
    def effective_minor_radius(self) -> float:
        """``a sqrt((1 + kappa**2) / 2)``, the radius of the equivalent circle.

        This is the radius at which a circular plasma would have the same
        cylindrical safety factor as the elongated one, and it is what the
        poloidal field average is taken over.
        """
        return self.minor_radius * math.sqrt((1.0 + self.elongation**2) / 2.0)
