"""Deuterium tritium fusion cross-section and Maxwellian reactivity.

Both the cross-section parameterisation and the thermal reactivity
parameterisation of Bosch and Hale (1992) are implemented. They are separate fits
in the same paper, derived from the same evaluated data by different procedures,
so evaluating the Maxwellian average of the cross-section fit numerically and
comparing it against the reactivity fit is an external check on this module that
does not depend on any number recalled from a table.

Reference:
    H.-S. Bosch and G. M. Hale, "Improved formulas for fusion cross-sections and
    thermal reactivities", Nuclear Fusion 32 (1992) 611, DOI 10.1088/0029-5515/32/4/I07.
    The reactivity coefficients are Table VII, the cross-section coefficients
    Table IV, and the Gamow constant and reduced mass energy Table I.
"""

from __future__ import annotations

import math
from typing import Final

import numpy as np
from scipy import integrate

from arc_benchmark._types import FloatArray
from arc_benchmark.model.constants import LIGHT_SPEED_CM_S

__all__ = [
    "REACTIVITY_TEMPERATURE_RANGE_KEV",
    "cross_section_barn",
    "dt_reactivity",
    "dt_reactivity_from_cross_section",
]

# Table I: Gamow constant in keV**(1/2) and reduced mass energy in keV, for T(d,n)4He.
_GAMOW_KEV_HALF: Final[float] = 34.3827
_REDUCED_MASS_KEV: Final[float] = 1124656.0

# Table VII, T(d,n)4He, valid for 0.2 keV to 100 keV with a stated maximum
# relative error of 0.25 percent against the underlying evaluation.
_C1: Final[float] = 1.17302e-9
_C2: Final[float] = 1.51361e-2
_C3: Final[float] = 7.51886e-2
_C4: Final[float] = 4.60643e-3
_C5: Final[float] = 1.35000e-2
_C6: Final[float] = -1.06750e-4
_C7: Final[float] = 1.36600e-5

REACTIVITY_TEMPERATURE_RANGE_KEV: Final[tuple[float, float]] = (0.2, 100.0)
"""The interval over which the Bosch and Hale reactivity fit is published."""

# Table IV, T(d,n)4He astrophysical S-function, valid for 0.5 keV to 550 keV.
_A: Final[tuple[float, float, float, float, float]] = (
    6.927e4,
    7.454e8,
    2.050e6,
    5.2002e4,
    0.0,
)
_B: Final[tuple[float, float, float, float]] = (6.38e1, -9.95e-1, 6.981e-5, 1.728e-4)

_CROSS_SECTION_RANGE_KEV: Final[tuple[float, float]] = (0.5, 550.0)
_MILLIBARN_PER_BARN: Final[float] = 1.0e-3
_SQUARE_CENTIMETRE_PER_BARN: Final[float] = 1.0e-24


def dt_reactivity(temperature_kev: FloatArray | float) -> FloatArray:
    """Maxwellian averaged deuterium tritium reactivity in cubic metre per second.

    Implements the Bosch and Hale (1992) closed form, equations (12) and (13),

        theta = T / (1 - T (C2 + T (C4 + T C6)) / (1 + T (C3 + T (C5 + T C7))))
        xi = (B_G**2 / (4 theta))**(1/3)
        <sigma v> = C1 theta sqrt(xi / (m_r c**2 T**3)) exp(-3 xi)

    with ``T`` in keV. The published fit returns cubic centimetre per second; the
    result is converted to cubic metre per second here because every other
    quantity in this package is SI.

    Args:
        temperature_kev: Ion temperature in keV. Scalars are accepted and
            returned as a zero-dimensional array.

    Returns:
        Reactivity in cubic metre per second, the same shape as the input.

    Raises:
        ValueError: If any temperature lies outside the published fit range.
    """
    t = np.asarray(temperature_kev, dtype=np.float64)
    low, high = REACTIVITY_TEMPERATURE_RANGE_KEV
    if np.any(t < low) or np.any(t > high):
        raise ValueError(
            f"temperature must lie in [{low}, {high}] keV where the Bosch and Hale fit is "
            f"published, got range [{np.min(t)}, {np.max(t)}]"
        )

    numerator = t * (_C2 + t * (_C4 + t * _C6))
    denominator = 1.0 + t * (_C3 + t * (_C5 + t * _C7))
    theta = t / (1.0 - numerator / denominator)
    xi = (_GAMOW_KEV_HALF**2 / (4.0 * theta)) ** (1.0 / 3.0)

    sigma_v_cm3 = _C1 * theta * np.sqrt(xi / (_REDUCED_MASS_KEV * t**3)) * np.exp(-3.0 * xi)
    return np.asarray(sigma_v_cm3 * 1.0e-6, dtype=np.float64)


def cross_section_barn(energy_kev: FloatArray | float) -> FloatArray:
    """Deuterium tritium fusion cross-section in barn against centre-of-mass energy.

    Implements the Bosch and Hale (1992) astrophysical S-function form,
    equations (8) and (9),

        S(E) = (A1 + E (A2 + E (A3 + E (A4 + E A5)))) / (1 + E (B1 + E (B2 + E (B3 + E B4))))
        sigma(E) = S(E) / (E exp(B_G / sqrt(E)))

    The published S-function is tabulated in keV millibarn, so the quotient above
    is a cross-section in millibarn. It is converted to barn on return, which is
    the unit the peak value is normally quoted in. The check that this factor is
    right is the location and height of the resonance: this function peaks at
    5.07 barn at a centre-of-mass energy of 64.8 keV, which is a deuteron
    laboratory energy of 108 keV, against the accepted 5 barn near 64 keV.

    Args:
        energy_kev: Centre-of-mass energy in keV.

    Returns:
        Cross-section in barn, the same shape as the input.

    Raises:
        ValueError: If any energy lies outside the published fit range.
    """
    e = np.asarray(energy_kev, dtype=np.float64)
    low, high = _CROSS_SECTION_RANGE_KEV
    if np.any(e < low) or np.any(e > high):
        raise ValueError(
            f"energy must lie in [{low}, {high}] keV where the Bosch and Hale cross-section "
            f"fit is published, got range [{np.min(e)}, {np.max(e)}]"
        )

    s_numerator = _A[0] + e * (_A[1] + e * (_A[2] + e * (_A[3] + e * _A[4])))
    s_denominator = 1.0 + e * (_B[0] + e * (_B[1] + e * (_B[2] + e * _B[3])))
    s_function = s_numerator / s_denominator
    millibarn = s_function / (e * np.exp(_GAMOW_KEV_HALF / np.sqrt(e)))
    return np.asarray(millibarn * _MILLIBARN_PER_BARN, dtype=np.float64)


def dt_reactivity_from_cross_section(
    temperature_kev: float,
    *,
    upper_energy_kev: float = 550.0,
    points: int = 4001,
) -> float:
    """Reactivity obtained by integrating the cross-section fit over a Maxwellian.

    This is the independent path to the same number that :func:`dt_reactivity`
    returns from the closed-form fit. It evaluates

        <sigma v> = c sqrt(8 / (pi m_r c**2)) T**(-3/2) integral sigma(E) E exp(-E / T) dE

    over centre-of-mass energy, with ``sigma`` in square centimetre and the
    result in cubic metre per second. The integrand is smooth and decays
    exponentially, so a composite Simpson rule on a fixed grid is used rather
    than an adaptive quadrature: a fixed grid gives the same answer on every
    machine, which an adaptive rule with a floating-point stopping test does not.

    Args:
        temperature_kev: Ion temperature in keV.
        upper_energy_kev: Upper limit of the energy integral, in keV. The default
            is the top of the published cross-section range.
        points: Number of grid points, forced odd for the Simpson rule.

    Returns:
        Reactivity in cubic metre per second.
    """
    if temperature_kev <= 0.0:
        raise ValueError(f"temperature must be positive, got {temperature_kev}")
    odd_points = points if points % 2 == 1 else points + 1
    low = _CROSS_SECTION_RANGE_KEV[0]
    energy = np.linspace(low, upper_energy_kev, odd_points, dtype=np.float64)

    sigma_cm2 = cross_section_barn(energy) * _SQUARE_CENTIMETRE_PER_BARN
    integrand = sigma_cm2 * energy * np.exp(-energy / temperature_kev)
    integral = float(integrate.simpson(integrand, x=energy))

    prefactor = LIGHT_SPEED_CM_S * math.sqrt(8.0 / (math.pi * _REDUCED_MASS_KEV))
    sigma_v_cm3 = float(prefactor * temperature_kev**-1.5 * integral)
    return sigma_v_cm3 * 1.0e-6
