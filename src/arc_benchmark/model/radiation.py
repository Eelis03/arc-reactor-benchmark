"""Radiated power densities: bremsstrahlung, synchrotron, and impurity line radiation.

Every function here returns a power density in watt per cubic metre and takes SI
densities in inverse cubic metre with temperatures in keV. None of them performs
any I/O or holds any state.

The three terms are not of equal quality. Bremsstrahlung is a closed-form result
with a small relativistic correction and is the most reliable. Synchrotron is a
Trubnikov asymptotic formula and is known to overestimate. Impurity line
radiation takes its atomic data from the caller because no atomic database ships
with this package.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np

from arc_benchmark._types import FloatArray
from arc_benchmark.model.constants import (
    ELECTRON_MASS,
    ELEMENTARY_CHARGE,
    JOULE_PER_KEV,
    LIGHT_SPEED,
    VACUUM_PERMITTIVITY,
)

__all__ = [
    "BREMSSTRAHLUNG_COEFFICIENT",
    "CORONAL_LOSS_PARAMETER",
    "LARMOR_COEFFICIENT",
    "TRUBNIKOV_OPACITY_COEFFICIENT",
    "ImpurityRadiator",
    "bremsstrahlung_density",
    "cyclotron_escape_factor",
    "line_radiation_density",
    "optically_thin_cyclotron_density",
    "synchrotron_density",
]

BREMSSTRAHLUNG_COEFFICIENT: Final[float] = 5.35e-37
"""Coefficient of the classical bremsstrahlung formula, SI with T in keV.

Wesson, "Tokamaks", gives ``P_br = 5.35e-37 Z_eff n_e**2 T_e**(1/2)`` in watt per
cubic metre for ``n_e`` in inverse cubic metre and ``T_e`` in keV.
"""

LARMOR_COEFFICIENT: Final[float] = (
    ELEMENTARY_CHARGE**4
    * JOULE_PER_KEV
    / (3.0 * math.pi * VACUUM_PERMITTIVITY * ELECTRON_MASS**3 * LIGHT_SPEED**3)
)
"""Coefficient of the optically thin cyclotron loss, assembled from constants.

Equals 1.5696e-14 in SI with the temperature in keV, so that
``P_thin = LARMOR_COEFFICIENT n_e B**2 T_keV`` in watt per cubic metre. Written
as an expression rather than as a literal because the literal usually quoted for
this, 6.2e-3 in megawatt per cubic metre with the density in units of 1e20, is
easy to transcribe with the wrong power of ten. The test suite checks that the
expression reproduces that literal to three significant figures, which is a check
on the arithmetic here and on the literal at the same time.
"""

TRUBNIKOV_OPACITY_COEFFICIENT: Final[float] = 6.04e3
"""Opacity coefficient of the Trubnikov escape factor, in metre keV per tesla.

This is the one number in the radiation model taken directly from a published
asymptotic result rather than derived here, and it is the least certain part of
the model. Its effect is confined to the synchrotron term.
"""

_ELECTRON_REST_ENERGY_KEV: Final[float] = 510.99895

CORONAL_LOSS_PARAMETER: Final[dict[str, float]] = {
    "beryllium": 0.0,
    "argon": 3.0e-33,
    "tungsten": 8.0e-32,
}
"""Representative coronal-equilibrium radiative loss parameters in watt cubic metre.

These are order-of-magnitude values for a core electron temperature between 5 and
20 keV, read from the published loss curves cited in the README. They exist so
that the examples have a defensible default, not as a substitute for an atomic
database. Beryllium is fully stripped above about 1 keV and therefore radiates
only bremsstrahlung, which the bremsstrahlung term already carries through
``Z_eff``; its line contribution is entered as zero for that reason. Any study
where impurity radiation matters should pass its own value.
"""


@dataclass(frozen=True, slots=True)
class ImpurityRadiator:
    """One impurity species, carrying both its charge and its atomic data.

    The charge and the concentration are what set fuel dilution and the effective
    charge in :mod:`arc_benchmark.algorithm.balance`, and the loss parameter is
    what sets line radiation here. They are held on one object so that a study
    cannot end up radiating from a species that the composition does not contain.

    Attributes:
        species: Name of the species, used only in reports.
        atomic_number: Nuclear charge ``Z``. The average charge state is taken as
            fully stripped, which holds for light species above about 1 keV and
            is an overestimate for tungsten.
        concentration: ``n_z / n_e``, dimensionless.
        loss_parameter: Radiative loss parameter ``L_z`` in watt cubic metre.
    """

    species: str
    atomic_number: int
    concentration: float
    loss_parameter: float

    def __post_init__(self) -> None:
        """Reject concentrations, charges, and loss parameters that are not physical."""
        if self.atomic_number < 1:
            raise ValueError(f"atomic_number must be at least 1, got {self.atomic_number}")
        if not 0.0 <= self.concentration <= 1.0:
            raise ValueError(f"concentration must lie in [0, 1], got {self.concentration}")
        if self.loss_parameter < 0.0:
            raise ValueError(f"loss_parameter must be non-negative, got {self.loss_parameter}")


def bremsstrahlung_density(
    electron_density: FloatArray | float,
    electron_temperature_kev: FloatArray | float,
    z_effective: float = 1.0,
    *,
    relativistic: bool = True,
) -> FloatArray:
    """Bremsstrahlung power density in watt per cubic metre.

    The classical result is

        P_br = 5.35e-37 Z_eff n_e**2 sqrt(T_e)

    with ``n_e`` in inverse cubic metre and ``T_e`` in keV. Above roughly 10 keV
    the electrons are relativistic enough for this to underestimate, and
    electron-electron bremsstrahlung, which vanishes non-relativistically, starts
    to contribute. Both are covered by the multiplier

        1 + 1.5 t + 0.7936 t + 1.874 t**2,  t = T_e / (m_e c**2)

    of Rider (1995), which raises the loss by 6 percent at 14 keV and by 13
    percent at 30 keV. Passing ``relativistic=False`` recovers the classical
    formula exactly, which is what the unit test checks against a hand
    evaluation.

    Args:
        electron_density: Electron density in inverse cubic metre.
        electron_temperature_kev: Electron temperature in keV.
        z_effective: Effective charge, dimensionless, at least one.
        relativistic: Whether to apply the relativistic and electron-electron
            correction.

    Returns:
        Power density in watt per cubic metre.

    Raises:
        ValueError: If ``z_effective`` is below one.
    """
    if z_effective < 1.0:
        raise ValueError(f"z_effective must be at least 1, got {z_effective}")

    n_e = np.asarray(electron_density, dtype=np.float64)
    t_e = np.asarray(electron_temperature_kev, dtype=np.float64)
    classical = BREMSSTRAHLUNG_COEFFICIENT * z_effective * n_e**2 * np.sqrt(t_e)
    if not relativistic:
        return np.asarray(classical, dtype=np.float64)

    t = t_e / _ELECTRON_REST_ENERGY_KEV
    correction = 1.0 + 1.5 * t + 0.7936 * t + 1.874 * t**2
    return np.asarray(classical * correction, dtype=np.float64)


def optically_thin_cyclotron_density(
    electron_density: FloatArray | float,
    electron_temperature_kev: FloatArray | float,
    toroidal_field: float,
) -> FloatArray:
    """Cyclotron emission with no reabsorption at all, in watt per cubic metre.

    The Larmor power of one electron gyrating in a field ``B`` is
    ``e**4 B**2 v_perp**2 / (6 pi eps_0 m**2 c**3)``. Averaging the perpendicular
    speed over a Maxwellian gives ``<v_perp**2> = 2 T / m`` and therefore

        P_thin = n_e e**4 B**2 T / (3 pi eps_0 m**3 c**3)

    with ``T`` in joule. This is not a loss the plasma actually suffers: a
    tokamak core is optically thick to the low harmonics and reabsorbs almost all
    of this. It is computed because it is exact, it fixes the scale, and it is
    what :func:`cyclotron_escape_factor` is a fraction of. At the ITER reference
    point it comes to 1.5 MW per cubic metre, about fifty times the synchrotron
    loss that actually escapes and about four times the fusion power density.

    Args:
        electron_density: Electron density in inverse cubic metre.
        electron_temperature_kev: Electron temperature in keV.
        toroidal_field: Vacuum toroidal field on axis in tesla.

    Returns:
        Power density in watt per cubic metre.
    """
    if toroidal_field <= 0.0:
        raise ValueError(f"toroidal_field must be positive, got {toroidal_field}")
    n_e = np.asarray(electron_density, dtype=np.float64)
    t_e = np.asarray(electron_temperature_kev, dtype=np.float64)
    return np.asarray(LARMOR_COEFFICIENT * n_e * toroidal_field**2 * t_e, dtype=np.float64)


def cyclotron_escape_factor(
    electron_density: FloatArray | float,
    electron_temperature_kev: FloatArray | float,
    toroidal_field: float,
    minor_radius: float,
    wall_reflectivity: float = 0.9,
) -> FloatArray:
    """Fraction of the cyclotron emission that escapes, Trubnikov asymptotic form.

    The plasma is optically thick to the low cyclotron harmonics and thin above a
    cutoff harmonic that rises with temperature and falls with density and size.
    Trubnikov's asymptotic evaluation of where that cutoff falls gives

        Phi = sqrt((1 - R_w) T_e B_0 / (6.04e3 a n_e20))

    with ``T_e`` in keV, ``B_0`` in tesla, ``a`` in metre, and ``n_e20`` in units
    of 1e20 inverse cubic metre. ``R_w`` is the fraction reflected by the first
    wall and returned to the plasma, where it is reabsorbed.

    The factor is around 0.02 for ITER parameters and 0.04 for ARC parameters,
    which is what turns a 1.5 MW per cubic metre emission into a loss of tens of
    kilowatt per cubic metre.

    Args:
        electron_density: Electron density in inverse cubic metre.
        electron_temperature_kev: Electron temperature in keV.
        toroidal_field: Vacuum toroidal field on axis in tesla.
        minor_radius: Plasma minor radius in metre.
        wall_reflectivity: Reflected fraction, in [0, 1).

    Returns:
        The escaping fraction, dimensionless.

    Raises:
        ValueError: If the reflectivity is outside [0, 1) or the geometry is not
            positive.
    """
    if not 0.0 <= wall_reflectivity < 1.0:
        raise ValueError(f"wall_reflectivity must lie in [0, 1), got {wall_reflectivity}")
    if toroidal_field <= 0.0:
        raise ValueError(f"toroidal_field must be positive, got {toroidal_field}")
    if minor_radius <= 0.0:
        raise ValueError(f"minor_radius must be positive, got {minor_radius}")

    n_e20 = np.asarray(electron_density, dtype=np.float64) / 1.0e20
    t_e = np.asarray(electron_temperature_kev, dtype=np.float64)
    opacity = TRUBNIKOV_OPACITY_COEFFICIENT * minor_radius * n_e20 / toroidal_field
    return np.asarray(np.sqrt((1.0 - wall_reflectivity) * t_e / opacity), dtype=np.float64)


def synchrotron_density(
    electron_density: FloatArray | float,
    electron_temperature_kev: FloatArray | float,
    toroidal_field: float,
    minor_radius: float,
    wall_reflectivity: float = 0.9,
) -> FloatArray:
    """Synchrotron power density in watt per cubic metre after reabsorption.

    The product of :func:`optically_thin_cyclotron_density` and
    :func:`cyclotron_escape_factor`, which collapses to

        P_syn proportional to n_e**0.5 T_e**1.5 B_0**2.5 a**(-0.5) sqrt(1 - R_w)

    The two and a half power of field is what makes this the price of the
    high-field pathway. Going from the 5.3 T of ITER to the 9.2 T of ARC
    multiplies this loss by a factor of four at fixed everything else, and going
    on to the 12.2 T of SPARC multiplies it by another factor of two.

    This is the weakest term in the model. Trubnikov's asymptotic form
    overestimates at high temperature, where the correct treatment needs the full
    emission profile rather than a single cutoff harmonic. The design notes
    record what the alternative would be and what it would cost.

    Args:
        electron_density: Electron density in inverse cubic metre.
        electron_temperature_kev: Electron temperature in keV.
        toroidal_field: Vacuum toroidal field on axis in tesla.
        minor_radius: Plasma minor radius in metre.
        wall_reflectivity: Reflected fraction, in [0, 1). A polished metal first
            wall is usually taken at 0.9.

    Returns:
        Power density in watt per cubic metre.
    """
    thin = optically_thin_cyclotron_density(
        electron_density, electron_temperature_kev, toroidal_field
    )
    escaping = cyclotron_escape_factor(
        electron_density,
        electron_temperature_kev,
        toroidal_field,
        minor_radius,
        wall_reflectivity,
    )
    return np.asarray(thin * escaping, dtype=np.float64)


def line_radiation_density(
    electron_density: FloatArray | float,
    radiators: tuple[ImpurityRadiator, ...],
) -> FloatArray:
    """Impurity line radiation power density in watt per cubic metre.

    Each species contributes ``n_e n_z L_z``, with ``n_z = c_z n_e``, so the sum
    is ``n_e**2 sum(c_z L_z)``. The loss parameter is supplied by the caller and
    is assumed to already be the coronal-equilibrium value at the core electron
    temperature; the weak temperature dependence of ``L_z`` for a high-Z species
    above 5 keV is not modelled.

    Args:
        electron_density: Electron density in inverse cubic metre.
        radiators: Impurity species present, possibly empty.

    Returns:
        Power density in watt per cubic metre.
    """
    n_e = np.asarray(electron_density, dtype=np.float64)
    total = sum(r.concentration * r.loss_parameter for r in radiators)
    return np.asarray(n_e**2 * total, dtype=np.float64)
