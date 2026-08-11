"""Parametric sweeps over field, size, and density.

A sweep is only meaningful once it says what it holds fixed, so the invariant is
part of the call signature and is recorded in the trace. Sweeping the toroidal
field while holding the density fixed and sweeping it while holding the beta
fixed are different questions with different answers, and the difference between
those answers is the whole argument for the high-field pathway.

Every sweep holds the cylindrical safety factor at its base value, so the plasma
current tracks the field and the size in the way an equilibrium requires. Without
that, a field sweep would be a sweep of the safety factor as well and the
measured exponents would mean nothing.
"""

from __future__ import annotations

import dataclasses
from collections.abc import Sequence
from enum import Enum

from arc_benchmark.algorithm.balance import LossPowerConvention, PlasmaState
from arc_benchmark.algorithm.constraints import ConstraintLimits, evaluate_constraints
from arc_benchmark.algorithm.operating import solve_operating_point
from arc_benchmark.algorithm.protocols import ConfinementScaling
from arc_benchmark.model.confinement import IPB98Y2
from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.model.limits import (
    current_at_fixed_q,
    cylindrical_safety_factor,
    greenwald_density,
)
from arc_benchmark.model.plant import PlantParameters, evaluate_plant
from arc_benchmark.pipeline.trace import SweepPoint, SweepTrace

__all__ = ["SweepInvariant", "density_sweep", "field_sweep", "radius_sweep"]


class SweepInvariant(Enum):
    """What the density is held to while the swept variable changes.

    ``FIXED_DENSITY`` holds the electron density at its base value. This isolates
    the effect of the swept variable on confinement alone, since fusion power at
    fixed density and temperature does not depend on the field at all.

    ``FIXED_GREENWALD_FRACTION`` holds the density at the same fraction of the
    Greenwald limit. Since the Greenwald limit tracks the current density and the
    current tracks the field at fixed safety factor, this makes the density
    proportional to the field.

    ``FIXED_BETA`` holds the toroidal beta, which at fixed temperature makes the
    density proportional to the square of the field. This is the branch on which
    the high-field argument is usually made.
    """

    FIXED_DENSITY = "fixed electron density"
    FIXED_GREENWALD_FRACTION = "fixed Greenwald fraction"
    FIXED_BETA = "fixed toroidal beta"


def _rebuilt_state(
    base: PlasmaState,
    geometry: PlasmaGeometry,
    toroidal_field: float,
    invariant: SweepInvariant,
    base_q: float,
    base_greenwald_fraction: float,
) -> PlasmaState:
    """Build the state at a new geometry and field, holding the safety factor."""
    current = current_at_fixed_q(geometry, toroidal_field, base_q)

    if invariant is SweepInvariant.FIXED_DENSITY:
        density = base.electron_density
    elif invariant is SweepInvariant.FIXED_GREENWALD_FRACTION:
        # The Greenwald fraction is defined on the line-averaged density, so the
        # target line average is converted back to the volume average this state
        # carries. The ratio is one for a flat density, which every sweep in this
        # package uses, and the division is here so that a peaked sweep would
        # hold the fraction the limit is actually written in.
        line_averaged = base_greenwald_fraction * greenwald_density(current, geometry.minor_radius)
        density = line_averaged / base.profile.line_average_ratio()
    else:
        # Toroidal beta is proportional to n T / B**2, and the temperature is not
        # swept, so holding beta means the density follows the square of the field.
        density = base.electron_density * (toroidal_field / base.toroidal_field) ** 2

    return dataclasses.replace(
        base,
        geometry=geometry,
        toroidal_field=toroidal_field,
        plasma_current_ma=current,
        electron_density=density,
    )


def _evaluate(
    value: float,
    state: PlasmaState,
    scaling: ConfinementScaling,
    convention: LossPowerConvention,
    limits: ConstraintLimits | None,
    plant: PlantParameters | None,
) -> SweepPoint:
    point = solve_operating_point(state, scaling, convention)
    report = evaluate_constraints(point, limits)
    plant_result = None
    if plant is not None:
        plant_result = evaluate_plant(
            fusion_power_mw=point.terms.fusion_power_mw,
            alpha_power_mw=point.terms.alpha_power_mw,
            neutron_power_mw=point.terms.neutron_power_mw,
            auxiliary_power_mw=max(point.auxiliary_power_mw, 0.0),
            parameters=plant,
        )
    return SweepPoint(value=value, point=point, constraints=report, plant=plant_result)


def field_sweep(
    base: PlasmaState,
    fields: Sequence[float],
    invariant: SweepInvariant = SweepInvariant.FIXED_GREENWALD_FRACTION,
    scaling: ConfinementScaling = IPB98Y2,
    convention: LossPowerConvention = LossPowerConvention.SEPARATRIX,
    limits: ConstraintLimits | None = None,
    plant: PlantParameters | None = None,
) -> SweepTrace:
    """Sweep the toroidal field on axis at fixed size, shape, and safety factor.

    Args:
        base: The reference state defining everything not swept.
        fields: Toroidal field values in tesla.
        invariant: What the density is held to. See :class:`SweepInvariant`.
        scaling: Confinement scaling.
        convention: Loss power convention.
        limits: Constraint thresholds.
        plant: Balance of plant parameters, or ``None`` to skip the electrical
            accounting.

    Returns:
        The sweep trace, in the order the fields were given.
    """
    base_q = cylindrical_safety_factor(base.geometry, base.toroidal_field, base.plasma_current_ma)
    base_fraction = base.line_averaged_density / greenwald_density(
        base.plasma_current_ma, base.geometry.minor_radius
    )
    points = tuple(
        _evaluate(
            value=field,
            state=_rebuilt_state(base, base.geometry, field, invariant, base_q, base_fraction),
            scaling=scaling,
            convention=convention,
            limits=limits,
            plant=plant,
        )
        for field in fields
    )
    return SweepTrace(
        variable="toroidal field on axis",
        units="T",
        policy=(
            f"{invariant.value}, fixed size and shape, cylindrical safety factor held at "
            f"{base_q:.3f}"
        ),
        scaling_name=scaling.name,
        points=points,
    )


def radius_sweep(
    base: PlasmaState,
    radii: Sequence[float],
    invariant: SweepInvariant = SweepInvariant.FIXED_GREENWALD_FRACTION,
    scaling: ConfinementScaling = IPB98Y2,
    convention: LossPowerConvention = LossPowerConvention.SEPARATRIX,
    limits: ConstraintLimits | None = None,
    plant: PlantParameters | None = None,
) -> SweepTrace:
    """Sweep the major radius at fixed aspect ratio, shape, field, and safety factor.

    The minor radius follows the major radius so that the aspect ratio is
    preserved, which means the volume goes as the cube of the radius and the
    plasma current goes as the radius.

    Args:
        base: The reference state.
        radii: Major radius values in metre.
        invariant: What the density is held to.
        scaling: Confinement scaling.
        convention: Loss power convention.
        limits: Constraint thresholds.
        plant: Balance of plant parameters, or ``None``.

    Returns:
        The sweep trace.
    """
    base_q = cylindrical_safety_factor(base.geometry, base.toroidal_field, base.plasma_current_ma)
    base_fraction = base.line_averaged_density / greenwald_density(
        base.plasma_current_ma, base.geometry.minor_radius
    )
    epsilon = base.geometry.inverse_aspect_ratio

    points = tuple(
        _evaluate(
            value=radius,
            state=_rebuilt_state(
                base,
                PlasmaGeometry(
                    major_radius=radius,
                    minor_radius=epsilon * radius,
                    elongation=base.geometry.elongation,
                    triangularity=base.geometry.triangularity,
                ),
                base.toroidal_field,
                invariant,
                base_q,
                base_fraction,
            ),
            scaling=scaling,
            convention=convention,
            limits=limits,
            plant=plant,
        )
        for radius in radii
    )
    return SweepTrace(
        variable="major radius",
        units="m",
        policy=(
            f"{invariant.value}, fixed aspect ratio {base.geometry.aspect_ratio:.3f}, fixed "
            f"field, cylindrical safety factor held at {base_q:.3f}"
        ),
        scaling_name=scaling.name,
        points=points,
    )


def density_sweep(
    base: PlasmaState,
    densities: Sequence[float],
    scaling: ConfinementScaling = IPB98Y2,
    convention: LossPowerConvention = LossPowerConvention.SEPARATRIX,
    limits: ConstraintLimits | None = None,
    plant: PlantParameters | None = None,
) -> SweepTrace:
    """Sweep the electron density with everything else held at its base value.

    Args:
        base: The reference state.
        densities: Electron densities in inverse cubic metre.
        scaling: Confinement scaling.
        convention: Loss power convention.
        limits: Constraint thresholds.
        plant: Balance of plant parameters, or ``None``.

    Returns:
        The sweep trace.
    """
    points = tuple(
        _evaluate(
            value=density,
            state=dataclasses.replace(base, electron_density=density),
            scaling=scaling,
            convention=convention,
            limits=limits,
            plant=plant,
        )
        for density in densities
    )
    return SweepTrace(
        variable="electron density",
        units="m^-3",
        policy="everything except density held at the base value",
        scaling_name=scaling.name,
        points=points,
    )
