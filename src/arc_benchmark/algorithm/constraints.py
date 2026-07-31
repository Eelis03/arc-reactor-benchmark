"""Operational limits turned into verdicts.

The model layer evaluates each published limit as a number. This module decides
whether a design point respects it, and reports the answer in a form that cannot
be mistaken for success: a point that violates a limit carries the violation, and
the aggregate report is not satisfied unless every check is.

The utilisation of a check is its distance to the boundary expressed so that one
is the boundary and above one is a violation, whichever direction the limit runs
in. That makes the checks comparable, which is what allows a sweep to say which
constraint binds first as a parameter is pushed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from enum import Enum

from arc_benchmark.algorithm.operating import OperatingPoint
from arc_benchmark.model.limits import (
    DEFAULT_BOOTSTRAP_COEFFICIENT,
    TROYON_COEFFICIENT,
    average_poloidal_field,
    bootstrap_fraction,
    cylindrical_safety_factor,
    greenwald_density,
    lh_threshold_power,
    normalised_beta,
    poloidal_beta,
    safety_factor_95,
    toroidal_beta,
)

__all__ = [
    "ConstraintCheck",
    "ConstraintLimits",
    "ConstraintReport",
    "ConstraintSense",
    "cylindrical_q",
    "evaluate_constraints",
]


class ConstraintSense(Enum):
    """Which side of the limit a quantity must stay on."""

    UPPER = "upper"
    LOWER = "lower"


@dataclass(frozen=True, slots=True)
class ConstraintCheck:
    """One limit, evaluated against one design point.

    Attributes:
        name: Short identifier, stable across releases so a trace can be joined
            on it.
        quantity: What was measured, in words.
        value: The measured value.
        limit: The limit it was measured against.
        sense: Whether the limit is an upper or a lower bound.
        units: Units of ``value`` and ``limit``.
        reference: Where the limit comes from.
    """

    name: str
    quantity: str
    value: float
    limit: float
    sense: ConstraintSense
    units: str
    reference: str

    @property
    def utilisation(self) -> float:
        """Fraction of the limit used. One is the boundary, above one violates.

        For an upper bound this is ``value / limit``; for a lower bound it is
        ``limit / value``, so that both read the same way.
        """
        if self.sense is ConstraintSense.UPPER:
            return math.inf if self.limit == 0.0 else self.value / self.limit
        return math.inf if self.value == 0.0 else self.limit / self.value

    @property
    def satisfied(self) -> bool:
        """True when the design point respects this limit."""
        return self.utilisation <= 1.0

    @property
    def margin(self) -> float:
        """How much room is left, as a fraction. Negative when violated."""
        return 1.0 - self.utilisation


@dataclass(frozen=True, slots=True)
class ConstraintLimits:
    """The thresholds a design point is judged against.

    Attributes:
        greenwald_fraction: Maximum ``n_e / n_G``. One is the empirical limit
            itself; conservative designs stay below 0.9.
        beta_n_limit: Maximum normalised beta in percent metre tesla per
            megaampere.
        q95_minimum: Minimum safety factor at the 95 percent flux surface. Two is
            the ideal kink boundary; operating designs target near three.
        bootstrap_fraction_maximum: Maximum fraction of the plasma current that
            may be bootstrap driven. One is the hard bound: a plasma cannot
            bootstrap more current than it carries.
        bootstrap_coefficient: Coefficient of the ``sqrt(eps) beta_p`` estimate.
        require_h_mode: Whether to require that the loss power exceeds the L to H
            transition threshold. Set this false when the operating point was
            solved with an L-mode scaling, where the threshold does not apply.
    """

    greenwald_fraction: float = 1.0
    beta_n_limit: float = TROYON_COEFFICIENT
    q95_minimum: float = 2.0
    bootstrap_fraction_maximum: float = 1.0
    bootstrap_coefficient: float = DEFAULT_BOOTSTRAP_COEFFICIENT
    require_h_mode: bool = True


@dataclass(frozen=True, slots=True)
class ConstraintReport:
    """Every check applied to one design point.

    Attributes:
        checks: The checks, in a fixed order so that reports line up.
    """

    checks: tuple[ConstraintCheck, ...]

    @property
    def satisfied(self) -> bool:
        """True only when every check is satisfied."""
        return all(check.satisfied for check in self.checks)

    @property
    def violations(self) -> tuple[ConstraintCheck, ...]:
        """The checks that are not satisfied, worst first."""
        failed = [check for check in self.checks if not check.satisfied]
        return tuple(sorted(failed, key=lambda c: -c.utilisation))

    @property
    def binding(self) -> ConstraintCheck:
        """The check closest to its boundary, whether or not it is violated.

        This is the constraint a sweep runs into first as a parameter is pushed
        in some direction, which is the question a design study actually asks.
        """
        return max(self.checks, key=lambda c: c.utilisation)

    def named(self, name: str) -> ConstraintCheck:
        """Return the check with this name.

        Raises:
            KeyError: If no check with that name is present.
        """
        for check in self.checks:
            if check.name == name:
                return check
        raise KeyError(f"no constraint named {name!r}; have {[c.name for c in self.checks]}")


def evaluate_constraints(
    point: OperatingPoint,
    limits: ConstraintLimits | None = None,
) -> ConstraintReport:
    """Evaluate every operational limit against a solved operating point.

    Args:
        point: A solved operating point.
        limits: Thresholds to judge against. Defaults to
            :class:`ConstraintLimits` with its published values.

    Returns:
        The full report. Nothing is filtered: satisfied checks are carried
        alongside violated ones so that the margin on each is visible.
    """
    thresholds = limits if limits is not None else ConstraintLimits()
    state = point.state
    geometry = state.geometry

    density_limit = greenwald_density(state.plasma_current_ma, geometry.minor_radius)
    beta_t = toroidal_beta(state.thermal_pressure_pa, state.toroidal_field)
    beta_n = normalised_beta(
        beta_t, geometry.minor_radius, state.toroidal_field, state.plasma_current_ma
    )
    q95 = safety_factor_95(geometry, state.toroidal_field, state.plasma_current_ma)
    b_poloidal = average_poloidal_field(geometry, state.plasma_current_ma)
    beta_p = poloidal_beta(state.thermal_pressure_pa, b_poloidal)
    f_bootstrap = bootstrap_fraction(geometry, beta_p, thresholds.bootstrap_coefficient)

    checks = [
        ConstraintCheck(
            name="greenwald",
            quantity="line-averaged electron density",
            value=state.electron_density,
            limit=thresholds.greenwald_fraction * density_limit,
            sense=ConstraintSense.UPPER,
            units="m^-3",
            reference="Greenwald et al., Nuclear Fusion 28 (1988) 2199",
        ),
        ConstraintCheck(
            name="troyon",
            quantity="normalised beta",
            value=beta_n,
            limit=thresholds.beta_n_limit,
            sense=ConstraintSense.UPPER,
            units="% m T / MA",
            reference="Troyon et al., Plasma Physics and Controlled Fusion 26 (1984) 209",
        ),
        ConstraintCheck(
            name="safety_factor",
            quantity="safety factor at the 95 percent flux surface",
            value=q95,
            limit=thresholds.q95_minimum,
            sense=ConstraintSense.LOWER,
            units="dimensionless",
            reference="ITER physics design guidelines, Uckan et al. (1990)",
        ),
        ConstraintCheck(
            name="bootstrap",
            quantity="bootstrap current fraction",
            value=f_bootstrap,
            limit=thresholds.bootstrap_fraction_maximum,
            sense=ConstraintSense.UPPER,
            units="dimensionless",
            reference="ITER Physics Basis Chapter 2, after Wilson, Nuclear Fusion 32 (1992) 257",
        ),
    ]

    if thresholds.require_h_mode:
        checks.append(
            ConstraintCheck(
                name="lh_threshold",
                quantity="power crossing the separatrix",
                value=point.loss_power_mw,
                limit=lh_threshold_power(
                    state.electron_density, state.toroidal_field, geometry.surface_area
                ),
                sense=ConstraintSense.LOWER,
                units="MW",
                reference="Martin, Takizuka et al., J. Phys. Conf. Ser. 123 (2008) 012033",
            )
        )

    return ConstraintReport(checks=tuple(checks))


def cylindrical_q(point: OperatingPoint) -> float:
    """Cylindrical safety factor of a solved point, for reporting alongside q95."""
    return cylindrical_safety_factor(
        point.state.geometry, point.state.toroidal_field, point.state.plasma_current_ma
    )
