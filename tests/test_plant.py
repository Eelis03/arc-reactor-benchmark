"""Tier one: the electrical accounting."""

from __future__ import annotations

import dataclasses

import pytest

from arc_benchmark.model.plant import PlantParameters, evaluate_plant

_PARAMETERS = PlantParameters()


def _evaluate(auxiliary_mw: float = 40.0, fusion_mw: float = 500.0) -> object:
    return evaluate_plant(
        fusion_power_mw=fusion_mw,
        alpha_power_mw=0.2 * fusion_mw,
        neutron_power_mw=0.8 * fusion_mw,
        auxiliary_power_mw=auxiliary_mw,
        parameters=_PARAMETERS,
    )


def test_thermal_power_is_the_stated_sum() -> None:
    """Thermal power is the multiplied neutron power plus everything the plasma absorbs.

    Tolerance: three multiplications and two additions, so 1e-14 relative.
    """
    result = evaluate_plant(500.0, 100.0, 400.0, 40.0, _PARAMETERS)
    expected = _PARAMETERS.blanket_multiplication * 400.0 + 100.0 + 40.0
    assert result.thermal_power_mw == pytest.approx(expected, rel=1.0e-14)


def test_gross_output_is_the_thermal_efficiency_times_thermal_power() -> None:
    """No hidden factor sits between the thermal power and the generator."""
    result = evaluate_plant(500.0, 100.0, 400.0, 40.0, _PARAMETERS)
    assert result.gross_electric_mw == pytest.approx(
        _PARAMETERS.thermal_efficiency * result.thermal_power_mw, rel=1.0e-14
    )


def test_recirculating_power_is_the_sum_of_its_named_parts() -> None:
    """Every draw is named and enters the total exactly once."""
    result = evaluate_plant(500.0, 100.0, 400.0, 40.0, _PARAMETERS)
    parts = (
        result.heating_draw_mw
        + result.cryoplant_mw
        + result.balance_of_plant_mw
        + result.tritium_and_auxiliary_mw
    )
    assert result.recirculating_mw == pytest.approx(parts, rel=1.0e-14)
    assert result.net_electric_mw == pytest.approx(
        result.gross_electric_mw - result.recirculating_mw, rel=1.0e-14
    )


def test_heating_draw_is_the_coupled_power_over_the_wall_plug_efficiency() -> None:
    """A wall plug efficiency of one half doubles the grid draw."""
    parameters = dataclasses.replace(_PARAMETERS, heating_wallplug_efficiency=0.5)
    result = evaluate_plant(500.0, 100.0, 400.0, 40.0, parameters)
    assert result.heating_draw_mw == pytest.approx(80.0, rel=1.0e-14)


def test_net_output_falls_as_the_heating_requirement_rises() -> None:
    """More auxiliary heating costs more grid power than the extra heat recovers.

    Auxiliary power enters the thermal balance once, weighted by the thermal
    efficiency of 0.40, and leaves the electrical balance once, weighted by the
    inverse of the wall plug efficiency of 0.55, which is 1.82. The second is
    larger, so net output falls.
    """
    low = evaluate_plant(500.0, 100.0, 400.0, 20.0, _PARAMETERS)
    high = evaluate_plant(500.0, 100.0, 400.0, 120.0, _PARAMETERS)
    assert high.net_electric_mw < low.net_electric_mw


def test_a_plant_can_consume_more_than_it_produces() -> None:
    """A large enough heating requirement drives the net output negative.

    This is the failure mode the whole efficiency study exists to find, so it is
    asserted to be reachable and reported rather than clipped at zero.
    """
    result = evaluate_plant(500.0, 100.0, 400.0, 400.0, _PARAMETERS)
    assert result.net_electric_mw < 0.0
    assert result.net_efficiency < 0.0
    assert result.recirculating_fraction > 1.0
    assert result.engineering_gain < 1.0


def test_net_efficiency_is_relative_to_fusion_power() -> None:
    """The reported efficiency divides net electricity by fusion power."""
    result = evaluate_plant(500.0, 100.0, 400.0, 40.0, _PARAMETERS)
    assert result.net_efficiency == pytest.approx(result.net_electric_mw / 500.0, rel=1.0e-14)


def test_blanket_multiplication_only_multiplies_the_neutron_power() -> None:
    """Alpha and auxiliary power reach the coolant unmultiplied.

    Multiplying them as well would be a double count: neither passes through the
    breeding blanket, so neither can benefit from the exothermic lithium capture.
    """
    plain = dataclasses.replace(_PARAMETERS, blanket_multiplication=1.0)
    boosted = dataclasses.replace(_PARAMETERS, blanket_multiplication=1.5)
    difference = (
        evaluate_plant(500.0, 100.0, 400.0, 40.0, boosted).thermal_power_mw
        - evaluate_plant(500.0, 100.0, 400.0, 40.0, plain).thermal_power_mw
    )
    assert difference == pytest.approx(0.5 * 400.0, rel=1.0e-14)


def test_plant_rejects_unphysical_parameters() -> None:
    """Efficiencies outside their ranges and negative loads all raise."""
    with pytest.raises(ValueError, match="thermal_efficiency"):
        PlantParameters(thermal_efficiency=1.2)
    with pytest.raises(ValueError, match="blanket_multiplication"):
        PlantParameters(blanket_multiplication=0.9)
    with pytest.raises(ValueError, match="heating_wallplug_efficiency"):
        PlantParameters(heating_wallplug_efficiency=1.5)
    with pytest.raises(ValueError, match="fixed electrical loads"):
        PlantParameters(cryoplant_mw=-1.0)
    with pytest.raises(ValueError, match="balance_of_plant_fraction"):
        PlantParameters(balance_of_plant_fraction=1.0)


def test_plant_rejects_negative_input_powers() -> None:
    """A negative power is a programming error, not an operating point."""
    with pytest.raises(ValueError, match="auxiliary_power_mw"):
        evaluate_plant(500.0, 100.0, 400.0, -1.0, _PARAMETERS)
    with pytest.raises(ValueError, match="fusion_power_mw"):
        evaluate_plant(-1.0, 100.0, 400.0, 40.0, _PARAMETERS)
