"""Lawson criterion, ignition condition, and the triple product.

The Lawson condition is the density-confinement product at which the plasma
sustains itself. It is derived here from the same power balance the rest of the
package solves, rather than quoted, so that the fuel dilution, the effective
charge, and the impurity content of a specific design point all enter it.

Only losses that scale as the square of density can appear in a
density-independent Lawson condition. Bremsstrahlung does, and so does impurity
line radiation in the coronal model used here. Synchrotron radiation does not: it
carries a half power of density, so including it would make the condition a
function of density and the curve would no longer be a curve. That is stated
rather than quietly handled, and the size of what is left out is reported by the
benchmark cases.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Final

import numpy as np
from scipy import optimize

from arc_benchmark.algorithm.balance import PlasmaComposition
from arc_benchmark.model.constants import ALPHA_ENERGY_KEV, ALPHA_FRACTION, JOULE_PER_KEV
from arc_benchmark.model.radiation import BREMSSTRAHLUNG_COEFFICIENT
from arc_benchmark.model.reactivity import REACTIVITY_TEMPERATURE_RANGE_KEV, dt_reactivity

__all__ = [
    "LawsonPoint",
    "lawson_n_tau",
    "lawson_triple_product",
    "optimum_lawson_temperature",
]

_ELECTRON_REST_ENERGY_KEV: Final[float] = 510.99895
_OPTIMUM_TOLERANCE_KEV: Final[float] = 1.0e-6


@dataclass(frozen=True, slots=True)
class LawsonPoint:
    """The Lawson requirement at one temperature.

    Attributes:
        temperature_kev: Temperature the condition was evaluated at.
        gain: Fusion gain the condition was evaluated for. Infinite for ignition.
        n_tau: Required ``n_e tau_E`` in inverse cubic metre second, infinite when
            radiation alone exceeds the alpha heating at this temperature.
        triple_product: Required ``n_e T tau_E`` in inverse cubic metre keV second.
    """

    temperature_kev: float
    gain: float
    n_tau: float
    triple_product: float


def lawson_n_tau(
    temperature_kev: float,
    composition: PlasmaComposition,
    gain: float = math.inf,
    *,
    include_line_radiation: bool = True,
    relativistic_bremsstrahlung: bool = True,
) -> float:
    """Required ``n_e tau_E`` in inverse cubic metre second.

    Setting sources equal to sinks with an auxiliary power of ``P_fus / Q``,

        g (f/2)**2 <sigma v> E_alpha n_e**2 = (3/2)(1 + n_i/n_e) n_e T / tau_E
                                              + (C_br Z_eff sqrt(T) + n_e**-0 L) n_e**2

    with ``g = 1 + 1 / (f_alpha Q)`` and ``f_alpha`` the alpha fraction of the
    reaction energy. Every term carries ``n_e**2`` except the transport term,
    which carries ``n_e``, so dividing through leaves a condition on ``n_e tau_E``
    alone.

    Args:
        temperature_kev: Temperature in keV, inside the Bosch and Hale range.
        composition: Fuel dilution, effective charge, and impurities.
        gain: Target fusion gain ``Q``. The default, infinity, is ignition.
        include_line_radiation: Whether the impurity line term enters the
            condition. It scales as the square of density, so it legitimately
            can.
        relativistic_bremsstrahlung: Whether to apply the relativistic and
            electron-electron correction to the bremsstrahlung term. Kept as a
            switch so that the classical Lawson curve can be reproduced exactly.

    Returns:
        Required ``n_e tau_E`` in inverse cubic metre second, or infinity when
        the radiation loss alone exceeds the heating at this temperature.

    Raises:
        ValueError: If the gain is not positive.
    """
    if gain <= 0.0:
        raise ValueError(f"gain must be positive, got {gain}")

    reactivity = float(dt_reactivity(temperature_kev))
    fuel_half = 0.5 * composition.fuel_fraction

    heating_multiplier = 1.0 + 1.0 / (ALPHA_FRACTION * gain) if math.isfinite(gain) else 1.0
    heating = (
        heating_multiplier * fuel_half**2 * reactivity * ALPHA_ENERGY_KEV * JOULE_PER_KEV
    )

    correction = 1.0
    if relativistic_bremsstrahlung:
        t = temperature_kev / _ELECTRON_REST_ENERGY_KEV
        correction = 1.0 + 1.5 * t + 0.7936 * t + 1.874 * t**2
    radiation = (
        BREMSSTRAHLUNG_COEFFICIENT
        * composition.z_effective
        * math.sqrt(temperature_kev)
        * correction
    )
    if include_line_radiation:
        radiation += sum(
            imp.concentration * imp.loss_parameter for imp in composition.impurities
        )

    net = heating - radiation
    if net <= 0.0:
        return math.inf

    stored = 1.5 * (1.0 + composition.ion_fraction) * temperature_kev * JOULE_PER_KEV
    return stored / net


def lawson_triple_product(
    temperature_kev: float,
    composition: PlasmaComposition,
    gain: float = math.inf,
    *,
    include_line_radiation: bool = True,
    relativistic_bremsstrahlung: bool = True,
) -> LawsonPoint:
    """Required ``n_e T tau_E`` at one temperature, with the inputs recorded."""
    n_tau = lawson_n_tau(
        temperature_kev,
        composition,
        gain,
        include_line_radiation=include_line_radiation,
        relativistic_bremsstrahlung=relativistic_bremsstrahlung,
    )
    return LawsonPoint(
        temperature_kev=temperature_kev,
        gain=gain,
        n_tau=n_tau,
        triple_product=n_tau * temperature_kev,
    )


def optimum_lawson_temperature(
    composition: PlasmaComposition,
    gain: float = math.inf,
    bracket_kev: tuple[float, float] = (5.0, 60.0),
    *,
    include_line_radiation: bool = True,
    relativistic_bremsstrahlung: bool = True,
) -> LawsonPoint:
    """Temperature minimising the required triple product, and that minimum.

    The triple product requirement falls with temperature while the reactivity
    rises faster than linearly, then rises again once the reactivity flattens
    above about 60 keV, so the requirement has an interior minimum. That minimum
    is the classical answer to what temperature a reactor should run at, before
    any engineering constraint is applied.

    The search uses a bounded Brent minimisation with an explicit absolute
    tolerance on the temperature. The tolerance is stated rather than defaulted
    because the returned temperature is pinned in the regression test and a
    tolerance that varies between SciPy releases would make that pin fragile.

    Args:
        composition: Fuel dilution, effective charge, and impurities.
        gain: Target fusion gain. The default, infinity, is ignition.
        bracket_kev: Interval to search, clipped to the Bosch and Hale range.
        include_line_radiation: Passed through to :func:`lawson_n_tau`.
        relativistic_bremsstrahlung: Passed through to :func:`lawson_n_tau`.

    Returns:
        The Lawson point at the minimising temperature.

    Raises:
        RuntimeError: If the minimiser reports that it did not converge.
    """
    low = max(bracket_kev[0], REACTIVITY_TEMPERATURE_RANGE_KEV[0])
    high = min(bracket_kev[1], REACTIVITY_TEMPERATURE_RANGE_KEV[1])
    if not low < high:
        raise ValueError(f"bracket does not intersect the reactivity range: {bracket_kev}")

    def objective(temperature_kev: float) -> float:
        value = lawson_n_tau(
            float(temperature_kev),
            composition,
            gain,
            include_line_radiation=include_line_radiation,
            relativistic_bremsstrahlung=relativistic_bremsstrahlung,
        )
        return value * float(temperature_kev) if math.isfinite(value) else np.inf

    result = optimize.minimize_scalar(
        objective,
        bounds=(low, high),
        method="bounded",
        options={"xatol": _OPTIMUM_TOLERANCE_KEV},
    )
    if not result.success:
        raise RuntimeError(f"triple product minimisation did not converge: {result.message}")

    return lawson_triple_product(
        float(result.x),
        composition,
        gain,
        include_line_radiation=include_line_radiation,
        relativistic_bremsstrahlung=relativistic_bremsstrahlung,
    )
