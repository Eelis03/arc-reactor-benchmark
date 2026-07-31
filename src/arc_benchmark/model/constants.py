"""Physical constants and deuterium tritium fusion energetics.

Constants are the 2018 CODATA values, which are exact for the elementary charge
and the speed of light under the 2019 revision of the SI. Reaction energies are
those of Bosch and Hale (1992), Table I, so that the energy release used here is
the one consistent with the reactivity fit in :mod:`arc_benchmark.model.reactivity`.
"""

from __future__ import annotations

from typing import Final

__all__ = [
    "ALPHA_ENERGY_KEV",
    "ALPHA_FRACTION",
    "DT_ENERGY_KEV",
    "ELECTRON_MASS",
    "ELEMENTARY_CHARGE",
    "JOULE_PER_KEV",
    "LIGHT_SPEED",
    "LIGHT_SPEED_CM_S",
    "NEUTRON_ENERGY_KEV",
    "NEUTRON_FRACTION",
    "VACUUM_PERMEABILITY",
    "VACUUM_PERMITTIVITY",
]

ELEMENTARY_CHARGE: Final[float] = 1.602176634e-19
"""Elementary charge in coulomb. Exact by the 2019 SI definition."""

JOULE_PER_KEV: Final[float] = 1.602176634e-16
"""Joule per kiloelectronvolt. Temperatures are carried in keV throughout."""

VACUUM_PERMEABILITY: Final[float] = 1.25663706212e-6
"""Magnetic constant in henry per metre. 2018 CODATA."""

VACUUM_PERMITTIVITY: Final[float] = 8.8541878128e-12
"""Electric constant in farad per metre. 2018 CODATA."""

ELECTRON_MASS: Final[float] = 9.1093837015e-31
"""Electron rest mass in kilogram. 2018 CODATA."""

LIGHT_SPEED: Final[float] = 2.99792458e8
"""Speed of light in metre per second. Exact by definition."""

LIGHT_SPEED_CM_S: Final[float] = LIGHT_SPEED * 100.0
"""Speed of light in centimetre per second.

Centimetre units appear only where the Bosch and Hale fit is evaluated, which is
published in cubic centimetre per second.
"""

DT_ENERGY_KEV: Final[float] = 17589.3
"""Energy released per deuterium tritium reaction in keV.

Bosch and Hale (1992) Table I give 17.589 MeV for T(d,n)4He.
"""

ALPHA_ENERGY_KEV: Final[float] = 3518.0
"""Kinetic energy of the alpha particle in keV, 3.518 MeV.

This is the fraction of the reaction energy that stays in the plasma and heats
it. The split follows from momentum conservation between a 4.0015 u alpha and a
1.0087 u neutron.
"""

NEUTRON_ENERGY_KEV: Final[float] = DT_ENERGY_KEV - ALPHA_ENERGY_KEV
"""Kinetic energy of the neutron in keV, 14.071 MeV. This leaves the plasma."""

ALPHA_FRACTION: Final[float] = ALPHA_ENERGY_KEV / DT_ENERGY_KEV
"""Fraction of fusion power deposited in the plasma, 0.2000."""

NEUTRON_FRACTION: Final[float] = NEUTRON_ENERGY_KEV / DT_ENERGY_KEV
"""Fraction of fusion power carried out by neutrons, 0.8000."""
