"""Steady-state operating point of the zero-dimensional balance.

The solve is closed form. In steady state the transport loss equals the stored
energy divided by the confinement time, and the confinement time is a power law
in that same loss power, so

    P_loss = W / (tau_1 P_loss**(-alpha))   ==>   P_loss = (W / tau_1)**(1 / (1 - alpha))

with ``tau_1`` the confinement time at one megawatt of loss power and ``alpha``
the power degradation exponent. There is no iteration, no tolerance, and no
convergence criterion, which is why the resulting operating points are safe to
pin in a regression test: the answer is a handful of floating-point operations
deep and is identical on any machine that implements IEEE 754.

The only iterative solve in this module is
:func:`solve_ignition_temperature`, which brackets a sign change and calls Brent's
method. That one converges to a stated tolerance and its convergence is asserted
before the result is returned.
"""

from __future__ import annotations

import dataclasses
import math
from dataclasses import dataclass
from typing import Final

from scipy import optimize

from arc_benchmark.algorithm.balance import (
    LossPowerConvention,
    PlasmaState,
    PowerTerms,
    power_terms,
)
from arc_benchmark.algorithm.protocols import ConfinementScaling

__all__ = [
    "OperatingPoint",
    "solve_ignition_temperature",
    "solve_operating_point",
]

_UNIT_LOSS_POWER_MW: Final[float] = 1.0
_IGNITION_TOLERANCE_KEV: Final[float] = 1.0e-9


@dataclass(frozen=True, slots=True)
class OperatingPoint:
    """A solved steady state, together with everything needed to audit it.

    Attributes:
        state: The plasma state that was solved.
        terms: Every state-dependent power term.
        scaling_name: Name of the confinement scaling used.
        convention: Which power the scaling was evaluated at.
        confinement_time_s: Energy confinement time in second, including the H
            factor.
        loss_power_mw: The power the confinement scaling was evaluated at.
        transport_power_mw: Stored energy divided by confinement time.
        auxiliary_power_mw: Heating power required to hold this state. Negative
            when alpha heating alone exceeds the losses, which is not a steady
            state at this temperature but is reported as the signed number rather
            than clipped, so that the balance closes exactly either way.
        fusion_gain: ``P_fus / P_aux``, infinite when the point is ignited.
        residual_mw: Sources minus sinks. Zero to rounding error by construction;
            it is carried so that a test can assert closure rather than trust it.
        triple_product: ``n_e T tau_E`` in inverse cubic metre keV second.
    """

    state: PlasmaState
    terms: PowerTerms
    scaling_name: str
    convention: LossPowerConvention
    confinement_time_s: float
    loss_power_mw: float
    transport_power_mw: float
    auxiliary_power_mw: float
    fusion_gain: float
    residual_mw: float
    triple_product: float

    @property
    def ignited(self) -> bool:
        """True when no external heating is needed to hold the balance."""
        return self.auxiliary_power_mw <= 0.0

    @property
    def radiated_fraction(self) -> float:
        """Radiated power divided by total heating power.

        Total heating is alpha plus auxiliary. Above one this would mean the
        plasma radiates more than it is heated by, which the solver never
        produces because the transport term is positive.
        """
        heating = self.terms.alpha_power_mw + self.auxiliary_power_mw
        if heating <= 0.0:
            return math.inf
        return self.terms.radiated_power_mw / heating


def solve_operating_point(
    state: PlasmaState,
    scaling: ConfinementScaling,
    convention: LossPowerConvention = LossPowerConvention.SEPARATRIX,
) -> OperatingPoint:
    """Solve the steady-state balance for the auxiliary power a state requires.

    Args:
        state: The plasma state, fully specified except for auxiliary heating.
        scaling: Any object satisfying
            :class:`~arc_benchmark.algorithm.protocols.ConfinementScaling`.
        convention: Which power the confinement scaling is evaluated at. See
            :class:`~arc_benchmark.algorithm.balance.LossPowerConvention`.

    Returns:
        The solved operating point.

    Raises:
        ValueError: If the scaling's power degradation is not strictly inside
            ``(0, 1)``, in which case no steady state exists.
    """
    alpha = scaling.power_degradation
    if not 0.0 < alpha < 1.0:
        raise ValueError(
            f"scaling {scaling.name!r} has power degradation {alpha}, which must lie strictly "
            "in (0, 1) for a steady state to exist"
        )

    terms = power_terms(state)
    tau_unit = state.confinement_multiplier * scaling.tau_e_at_unit_power(
        state.confinement_inputs(_UNIT_LOSS_POWER_MW)
    )

    loss_power_mw = (terms.stored_energy_mj / tau_unit) ** (1.0 / (1.0 - alpha))
    confinement_time_s = tau_unit * loss_power_mw**-alpha
    transport_power_mw = terms.stored_energy_mj / confinement_time_s

    radiated_counted = (
        terms.radiated_power_mw if convention is LossPowerConvention.SEPARATRIX else 0.0
    )
    auxiliary_power_mw = transport_power_mw + radiated_counted - terms.alpha_power_mw

    sources = terms.alpha_power_mw + auxiliary_power_mw
    sinks = radiated_counted + transport_power_mw
    residual_mw = sources - sinks

    fusion_gain = (
        terms.fusion_power_mw / auxiliary_power_mw if auxiliary_power_mw > 0.0 else math.inf
    )

    return OperatingPoint(
        state=state,
        terms=terms,
        scaling_name=scaling.name,
        convention=convention,
        confinement_time_s=confinement_time_s,
        loss_power_mw=loss_power_mw,
        transport_power_mw=transport_power_mw,
        auxiliary_power_mw=auxiliary_power_mw,
        fusion_gain=fusion_gain,
        residual_mw=residual_mw,
        triple_product=state.electron_density * state.temperature_kev * confinement_time_s,
    )


def solve_ignition_temperature(
    state: PlasmaState,
    scaling: ConfinementScaling,
    convention: LossPowerConvention = LossPowerConvention.SEPARATRIX,
    bracket_kev: tuple[float, float] = (2.0, 60.0),
) -> float | None:
    """Temperature at which the required auxiliary power crosses zero.

    Everything except temperature is held at the values carried by ``state``.
    Below the returned temperature the plasma needs external heating; above it
    the alpha power exceeds the losses.

    Brent's method is used on the bracket, after checking that the required
    auxiliary power actually changes sign across it. A bracket with no sign
    change means this plasma does not ignite anywhere in the interval, which is a
    result and is returned as ``None`` rather than as a number from a solve that
    did not converge.

    Args:
        state: Plasma state whose temperature is to be varied.
        scaling: Confinement scaling.
        convention: Loss power convention.
        bracket_kev: Temperature interval to search, in keV.

    Returns:
        The ignition temperature in keV, or ``None`` if the bracket contains no
        sign change.

    Raises:
        RuntimeError: If the root finder reports that it did not converge.
    """
    low, high = bracket_kev
    if not 0.0 < low < high:
        raise ValueError(f"bracket must satisfy 0 < low < high, got {bracket_kev}")

    def required_auxiliary(temperature_kev: float) -> float:
        trial = dataclasses.replace(state, temperature_kev=temperature_kev)
        return solve_operating_point(trial, scaling, convention).auxiliary_power_mw

    at_low = required_auxiliary(low)
    at_high = required_auxiliary(high)
    if at_low * at_high > 0.0:
        return None

    root, result = optimize.brentq(
        required_auxiliary,
        low,
        high,
        xtol=_IGNITION_TOLERANCE_KEV,
        full_output=True,
    )
    if not result.converged:
        raise RuntimeError(
            f"ignition temperature solve did not converge on bracket {bracket_kev}"
        )
    return float(root)
