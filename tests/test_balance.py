"""Tier one: the power balance and the operating point solve.

Closure, monotonicity, the ignition condition, and the field scaling exponents
are all properties the solver has to satisfy by construction, so each is asserted
directly rather than inferred from an output that happens to look right.
"""

from __future__ import annotations

import dataclasses
import itertools
import math

import pytest
from scipy import optimize

from arc_benchmark.algorithm.balance import (
    LossPowerConvention,
    PlasmaComposition,
    PlasmaState,
    power_terms,
)
from arc_benchmark.algorithm.operating import solve_ignition_temperature, solve_operating_point
from arc_benchmark.model.confinement import CONFINEMENT_SCALINGS, IPB98Y2
from arc_benchmark.model.constants import ALPHA_FRACTION, NEUTRON_FRACTION
from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.model.limits import current_at_fixed_q, cylindrical_safety_factor
from arc_benchmark.model.profiles import ProfileShape
from arc_benchmark.model.radiation import ImpurityRadiator
from arc_benchmark.pipeline.machines import machine

_ARC = machine("ARC").state


def _with(**changes: object) -> PlasmaState:
    return dataclasses.replace(_ARC, **changes)


@pytest.mark.parametrize("scaling_name", sorted(CONFINEMENT_SCALINGS))
@pytest.mark.parametrize("convention", list(LossPowerConvention))
def test_the_power_balance_closes(
    scaling_name: str, convention: LossPowerConvention
) -> None:
    """Sources equal sinks at the solved operating point.

    Tolerance: the residual is a difference of terms of order 100 MW, formed from
    about ten floating-point operations, so it cannot exceed roughly
    ``10 * 100 MW * 2.2e-16``, which is 2e-13 MW. The absolute tolerance is set
    to 1e-9 MW, four orders above that bound and still fifteen orders below the
    smallest term in the balance. It is derived from the magnitude of the terms
    and the machine epsilon, not from the residual that happens to be observed,
    which is exactly zero.
    """
    point = solve_operating_point(_ARC, CONFINEMENT_SCALINGS[scaling_name], convention)
    radiated = (
        point.terms.radiated_power_mw
        if convention is LossPowerConvention.SEPARATRIX
        else 0.0
    )
    sources = point.terms.alpha_power_mw + point.auxiliary_power_mw
    sinks = radiated + point.transport_power_mw
    assert sources - sinks == pytest.approx(0.0, abs=1.0e-9)
    assert point.residual_mw == pytest.approx(0.0, abs=1.0e-9)


def test_the_confinement_time_is_the_scaling_evaluated_at_the_solved_loss_power() -> None:
    """The solved point sits on the scaling, not near it.

    The closed-form inversion is only valid if the confinement time it produces
    is the one the scaling returns at the loss power it produces, so that
    identity is checked rather than trusted.
    """
    point = solve_operating_point(_ARC, IPB98Y2)
    predicted = _ARC.confinement_multiplier * IPB98Y2.tau_e(
        _ARC.confinement_inputs(point.loss_power_mw)
    )
    assert point.confinement_time_s == pytest.approx(predicted, rel=1.0e-13)
    assert point.transport_power_mw == pytest.approx(point.loss_power_mw, rel=1.0e-13)


def test_a_flat_density_makes_the_two_averages_identical() -> None:
    """The zero-dimensional default hands the scaling the density it carries.

    This is the property that keeps the line-average correction from touching any
    flat-profile result: the ratio is exactly one, so the correction is a
    multiplication by one and not an approximation applied everywhere.
    """
    assert _ARC.profile.is_flat
    assert _ARC.line_averaged_density == _ARC.electron_density
    assert _ARC.confinement_inputs(100.0).line_averaged_density_e19 == (
        _ARC.electron_density / 1.0e19
    )


def test_a_peaked_density_reaches_the_scaling_as_a_line_average() -> None:
    """The scaling is fed the chord average, which is what it was fitted against.

    The confinement scalings are regressions against an interferometer
    measurement, and a volume-averaged model owes them the conversion. The check
    is made twice: once on the input handed over, and once on the confinement
    time itself, which must move by the density exponent of the scaling and by
    nothing else, since only the density changed.
    """
    peaked = _with(profile=ProfileShape(0.4, 1.0))
    ratio = peaked.profile.line_average_ratio()
    assert ratio > 1.0

    assert peaked.line_averaged_density == pytest.approx(
        ratio * peaked.electron_density, rel=1.0e-15
    )
    assert peaked.confinement_inputs(100.0).line_averaged_density_e19 == pytest.approx(
        ratio * peaked.electron_density / 1.0e19, rel=1.0e-15
    )

    flat_tau = IPB98Y2.tau_e(_ARC.confinement_inputs(100.0))
    peaked_tau = IPB98Y2.tau_e(peaked.confinement_inputs(100.0))
    assert peaked_tau / flat_tau == pytest.approx(
        ratio**IPB98Y2.exponents.density, rel=1.0e-13
    )


def test_alpha_and_neutron_powers_partition_the_fusion_power() -> None:
    """The two products carry the whole reaction energy and nothing else."""
    terms = power_terms(_ARC)
    assert terms.alpha_power_mw + terms.neutron_power_mw == pytest.approx(
        terms.fusion_power_mw, rel=1.0e-14
    )
    assert pytest.approx(1.0, rel=1.0e-15) == ALPHA_FRACTION + NEUTRON_FRACTION
    assert terms.alpha_power_mw / terms.fusion_power_mw == pytest.approx(
        ALPHA_FRACTION, rel=1.0e-14
    )


def test_gain_increases_monotonically_with_confinement_time() -> None:
    """Longer confinement means less heating and therefore a higher gain.

    Everything except the H factor is held fixed, so the confinement time is the
    only thing that changes. Both the confinement time and the gain must rise
    without exception.
    """
    multipliers = [1.2 + 0.1 * step for step in range(0, 15)]
    solved = [
        solve_operating_point(_with(confinement_multiplier=m), IPB98Y2) for m in multipliers
    ]
    times = [p.confinement_time_s for p in solved]
    gains = [p.fusion_gain for p in solved]
    assert all(b > a for a, b in itertools.pairwise(times))
    assert all(b > a for a, b in itertools.pairwise(gains))


def test_ignition_is_exactly_where_alpha_heating_equals_the_losses() -> None:
    """At the ignition boundary the alpha power equals radiation plus transport.

    The boundary is located by root finding on the H factor, which is the only
    free parameter that moves the required auxiliary power monotonically. The
    identity that follows is a construction rather than a coincidence: the
    balance closes for every state, so an auxiliary power of zero forces the
    alpha power to carry every sink on its own.

    Tolerance: Brent's method is asked for a root to 1e-12 in the H factor, and
    the required auxiliary power changes by about 100 MW per unit of H near the
    root, so the residual auxiliary power at the returned root is of order
    1e-10 MW. The absolute tolerance is set at 1e-6 MW, four orders above that.
    """

    def required_auxiliary(multiplier: float) -> float:
        return solve_operating_point(
            _with(confinement_multiplier=multiplier), IPB98Y2
        ).auxiliary_power_mw

    root = optimize.brentq(required_auxiliary, 2.0, 8.0, xtol=1.0e-12)
    point = solve_operating_point(_with(confinement_multiplier=float(root)), IPB98Y2)

    assert point.auxiliary_power_mw == pytest.approx(0.0, abs=1.0e-6)
    assert point.terms.alpha_power_mw == pytest.approx(
        point.terms.radiated_power_mw + point.transport_power_mw, abs=1.0e-6
    )
    assert point.ignited
    assert math.isinf(point.fusion_gain)


def test_beyond_ignition_the_gain_is_reported_as_infinite() -> None:
    """A state past ignition returns an infinite gain and a negative auxiliary power."""
    point = solve_operating_point(_with(confinement_multiplier=8.0), IPB98Y2)
    assert point.auxiliary_power_mw < 0.0
    assert point.ignited
    assert math.isinf(point.fusion_gain)


def test_ignition_temperature_search_reports_no_root_rather_than_inventing_one() -> None:
    """A bracket whose endpoints share a sign returns ``None``.

    Two separate facts are asserted here. The first is that the ARC state as
    published does not ignite anywhere in the default bracket: with the
    IPB98(y,2) power degradation the transport loss rises as the temperature to
    the power 3.23 while the reactivity flattens above 30 keV, so raising the
    temperature never closes the balance on alpha heating alone.

    The second is that even a state which does ignite has a bounded ignition
    window rather than a half line. At an H factor of 8 the ARC state needs
    external heating below about 10 keV and above about 30 keV, and ignites in
    between. The default bracket therefore has the same sign at both ends and is
    correctly reported as containing no bracketed root, while a bracket placed
    across the lower crossing returns it.
    """
    assert solve_ignition_temperature(_ARC, IPB98Y2) is None

    ignitable = _with(confinement_multiplier=8.0)

    def required_auxiliary(temperature_kev: float) -> float:
        return solve_operating_point(
            dataclasses.replace(ignitable, temperature_kev=temperature_kev), IPB98Y2
        ).auxiliary_power_mw

    assert required_auxiliary(4.0) > 0.0
    assert required_auxiliary(14.0) < 0.0
    assert required_auxiliary(50.0) > 0.0
    assert solve_ignition_temperature(ignitable, IPB98Y2, bracket_kev=(4.0, 50.0)) is None

    root = solve_ignition_temperature(ignitable, IPB98Y2, bracket_kev=(4.0, 14.0))
    assert root is not None
    assert 4.0 < root < 14.0
    assert required_auxiliary(root) == pytest.approx(0.0, abs=1.0e-6)


def test_fusion_power_scales_as_the_fourth_power_of_field_at_fixed_beta() -> None:
    """Doubling the field at fixed beta multiplies fusion power by sixteen.

    Holding the toroidal beta at a fixed temperature makes the density go as the
    square of the field, and fusion power goes as the square of the density, so
    the analytic expectation is exactly ``2**4``. The safety factor is held by
    scaling the current with the field, so the equilibrium does not drift.

    Tolerance: the computation is a handful of multiplications and one call to
    the reactivity fit, which is identical at both points because the temperature
    is unchanged. 1e-12 relative is four orders above the accumulated rounding.
    """
    base_q = cylindrical_safety_factor(_ARC.geometry, _ARC.toroidal_field, _ARC.plasma_current_ma)
    doubled_field = 2.0 * _ARC.toroidal_field
    doubled = _with(
        toroidal_field=doubled_field,
        plasma_current_ma=current_at_fixed_q(_ARC.geometry, doubled_field, base_q),
        electron_density=4.0 * _ARC.electron_density,
    )
    ratio = power_terms(doubled).fusion_power_mw / power_terms(_ARC).fusion_power_mw
    assert ratio == pytest.approx(2.0**4, rel=1.0e-12)


def test_fusion_power_scales_as_the_square_of_field_at_fixed_greenwald_fraction() -> None:
    """Doubling the field at a fixed Greenwald fraction multiplies fusion power by four.

    The Greenwald limit tracks the current density and the current tracks the
    field at fixed safety factor, so the density goes as the field and fusion
    power as its square.
    """
    base_q = cylindrical_safety_factor(_ARC.geometry, _ARC.toroidal_field, _ARC.plasma_current_ma)
    doubled_field = 2.0 * _ARC.toroidal_field
    doubled = _with(
        toroidal_field=doubled_field,
        plasma_current_ma=current_at_fixed_q(_ARC.geometry, doubled_field, base_q),
        electron_density=2.0 * _ARC.electron_density,
    )
    ratio = power_terms(doubled).fusion_power_mw / power_terms(_ARC).fusion_power_mw
    assert ratio == pytest.approx(2.0**2, rel=1.0e-12)


def test_fusion_power_does_not_depend_on_field_at_fixed_density() -> None:
    """At fixed density, temperature, and volume, the field does not enter at all."""
    ratio = power_terms(_with(toroidal_field=2.0 * _ARC.toroidal_field)).fusion_power_mw
    assert ratio == pytest.approx(power_terms(_ARC).fusion_power_mw, rel=1.0e-14)


def test_fusion_power_scales_as_the_square_of_density_and_with_volume() -> None:
    """Unit consistency: the density and volume dependences are exactly as written."""
    base = power_terms(_ARC).fusion_power_mw
    denser = power_terms(_with(electron_density=2.0 * _ARC.electron_density)).fusion_power_mw
    assert denser == pytest.approx(4.0 * base, rel=1.0e-12)

    bigger_geometry = PlasmaGeometry(
        major_radius=2.0 * _ARC.geometry.major_radius,
        minor_radius=_ARC.geometry.minor_radius,
        elongation=_ARC.geometry.elongation,
        triangularity=_ARC.geometry.triangularity,
    )
    bigger = power_terms(_with(geometry=bigger_geometry)).fusion_power_mw
    assert bigger == pytest.approx(2.0 * base, rel=1.0e-12)


def test_stored_energy_scales_with_density_temperature_and_volume() -> None:
    """Unit consistency of the stored energy, one variable at a time."""
    base = _ARC.stored_energy_mj
    assert _with(electron_density=3.0 * _ARC.electron_density).stored_energy_mj == pytest.approx(
        3.0 * base, rel=1.0e-13
    )
    assert _with(temperature_kev=2.0 * _ARC.temperature_kev).stored_energy_mj == pytest.approx(
        2.0 * base, rel=1.0e-13
    )


def test_composition_enforces_quasineutrality_and_the_effective_charge() -> None:
    """Fuel dilution and effective charge follow from the impurity content.

    A pure hydrogenic plasma has a fuel fraction of one and an effective charge
    of one. Adding helium at a fraction ``f`` removes ``2 f`` of fuel and adds
    ``2 f`` to the effective charge. Adding a species of charge ``Z`` at
    concentration ``c`` removes ``Z c`` of fuel and adds ``c Z (Z - 1)``.
    """
    clean = PlasmaComposition()
    assert clean.fuel_fraction == pytest.approx(1.0, rel=1.0e-15)
    assert clean.z_effective == pytest.approx(1.0, rel=1.0e-15)

    ashy = PlasmaComposition(helium_fraction=0.05)
    assert ashy.fuel_fraction == pytest.approx(0.90, rel=1.0e-14)
    assert ashy.z_effective == pytest.approx(1.10, rel=1.0e-14)

    dirty = PlasmaComposition(
        helium_fraction=0.04,
        impurities=(ImpurityRadiator("beryllium", 4, 0.02, 0.0),),
    )
    assert dirty.fuel_fraction == pytest.approx(1.0 - 0.08 - 0.08, rel=1.0e-14)
    assert dirty.z_effective == pytest.approx(1.0 + 0.08 + 0.02 * 4 * 3, rel=1.0e-14)


def test_composition_rejects_a_plasma_with_no_fuel_left() -> None:
    """An impurity load that consumes the whole electron density raises."""
    with pytest.raises(ValueError, match="no fuel"):
        PlasmaComposition(impurities=(ImpurityRadiator("tungsten", 74, 0.02, 0.0),))


def test_dilution_reduces_fusion_power_by_the_square_of_the_fuel_fraction() -> None:
    """Fusion power carries the square of the fuel fraction, not the first power."""
    clean = _with(composition=PlasmaComposition())
    ashy = _with(composition=PlasmaComposition(helium_fraction=0.05))
    ratio = power_terms(ashy).fusion_power_mw / power_terms(clean).fusion_power_mw
    assert ratio == pytest.approx(0.90**2, rel=1.0e-13)


def test_the_two_loss_power_conventions_differ_by_the_radiated_power() -> None:
    """Subtracting radiation moves the required auxiliary power by exactly that much.

    The loss power itself is the same under both conventions, because it is fixed
    by the stored energy and the scaling alone. What changes is how much of the
    heating is accounted to it.
    """
    separatrix = solve_operating_point(_ARC, IPB98Y2, LossPowerConvention.SEPARATRIX)
    total = solve_operating_point(_ARC, IPB98Y2, LossPowerConvention.TOTAL)
    assert separatrix.loss_power_mw == pytest.approx(total.loss_power_mw, rel=1.0e-14)
    difference = separatrix.auxiliary_power_mw - total.auxiliary_power_mw
    assert difference == pytest.approx(separatrix.terms.radiated_power_mw, rel=1.0e-12)


def test_plasma_state_rejects_unphysical_values() -> None:
    """Non-positive density, temperature, field, current, or H factor all raise."""
    for name in (
        "toroidal_field",
        "plasma_current_ma",
        "electron_density",
        "temperature_kev",
        "confinement_multiplier",
    ):
        with pytest.raises(ValueError, match=name):
            _with(**{name: 0.0})
