"""Balance of plant: thermal conversion, recirculating power, and net efficiency.

This is the layer where a fusion power number becomes an electricity number, and
it is where a plasma that looks impressive can stop being a power plant. All
powers are in megawatt.

The accounting is deliberately explicit. Every term that consumes electricity is
named and enters the total once, so that a sweep can attribute a change in net
efficiency to a specific term rather than to an opaque recirculating fraction.
"""

from __future__ import annotations

from dataclasses import dataclass

__all__ = ["PlantParameters", "PlantResult", "evaluate_plant"]


@dataclass(frozen=True, slots=True)
class PlantParameters:
    """Engineering parameters of the power conversion and auxiliary systems.

    Attributes:
        thermal_efficiency: Gross electrical output divided by thermal power.
            0.40 is a steam cycle at fusion blanket outlet temperatures; a helium
            or molten salt Brayton cycle reaches 0.45 or above.
        blanket_multiplication: Thermal energy deposited in the blanket per unit
            of neutron energy entering it. Above one because lithium 6 neutron
            capture is exothermic at 4.78 MeV per reaction and beryllium or lead
            multipliers add further neutrons. A tritium breeding ratio near 1.1
            in a lithium blanket gives roughly 1.3.
        heating_wallplug_efficiency: Power coupled to the plasma divided by grid
            power drawn by the heating and current drive system. Radio frequency
            systems are usually quoted between 0.4 and 0.6, neutral beams lower.
        cryoplant_mw: Electrical power for the magnet cryoplant. A rare earth
            barium copper oxide magnet operating at 20 K instead of 4 K is the
            reason an ARC-class device can enter a small number here.
        balance_of_plant_fraction: Pumps, controls, and auxiliaries, as a
            fraction of gross electrical output.
        tritium_and_auxiliary_mw: Tritium plant, vacuum, and site load in
            megawatt.
    """

    thermal_efficiency: float = 0.40
    blanket_multiplication: float = 1.30
    heating_wallplug_efficiency: float = 0.55
    cryoplant_mw: float = 10.0
    balance_of_plant_fraction: float = 0.04
    tritium_and_auxiliary_mw: float = 5.0

    def __post_init__(self) -> None:
        """Reject parameter values that are not physically meaningful."""
        if not 0.0 < self.thermal_efficiency < 1.0:
            raise ValueError(
                f"thermal_efficiency must lie in (0, 1), got {self.thermal_efficiency}"
            )
        if self.blanket_multiplication < 1.0:
            raise ValueError(
                "blanket_multiplication must be at least 1, got "
                f"{self.blanket_multiplication}"
            )
        if not 0.0 < self.heating_wallplug_efficiency <= 1.0:
            raise ValueError(
                "heating_wallplug_efficiency must lie in (0, 1], got "
                f"{self.heating_wallplug_efficiency}"
            )
        if self.cryoplant_mw < 0.0 or self.tritium_and_auxiliary_mw < 0.0:
            raise ValueError("fixed electrical loads must be non-negative")
        if not 0.0 <= self.balance_of_plant_fraction < 1.0:
            raise ValueError(
                "balance_of_plant_fraction must lie in [0, 1), got "
                f"{self.balance_of_plant_fraction}"
            )


@dataclass(frozen=True, slots=True)
class PlantResult:
    """The electrical accounting of one operating point, all powers in megawatt.

    Attributes:
        thermal_power_mw: Heat delivered to the power conversion system.
        gross_electric_mw: Electrical output at the generator terminals.
        heating_draw_mw: Grid power drawn by heating and current drive.
        cryoplant_mw: Grid power drawn by the cryoplant.
        balance_of_plant_mw: Grid power drawn by pumps and auxiliaries.
        tritium_and_auxiliary_mw: Grid power drawn by the tritium plant and site.
        recirculating_mw: Sum of every electrical draw.
        net_electric_mw: Gross output less recirculating power. Negative when the
            plant consumes more electricity than it produces.
        recirculating_fraction: Recirculating power divided by gross output.
        net_efficiency: Net electrical output divided by fusion power.
        engineering_gain: Gross electrical output divided by recirculating power,
            the quantity usually written ``Q_eng``.
    """

    thermal_power_mw: float
    gross_electric_mw: float
    heating_draw_mw: float
    cryoplant_mw: float
    balance_of_plant_mw: float
    tritium_and_auxiliary_mw: float
    recirculating_mw: float
    net_electric_mw: float
    recirculating_fraction: float
    net_efficiency: float
    engineering_gain: float


def evaluate_plant(
    fusion_power_mw: float,
    alpha_power_mw: float,
    neutron_power_mw: float,
    auxiliary_power_mw: float,
    parameters: PlantParameters,
) -> PlantResult:
    """Convert a plasma power balance into an electrical balance.

    The thermal power is

        P_th = M P_n + P_alpha + P_aux

    because in steady state everything the plasma absorbs leaves it again as
    radiation or transported heat and is caught by the first wall and divertor,
    so alpha heating and auxiliary heating both reach the coolant. Only the
    neutron power is multiplied, since only it passes through the breeding
    blanket.

    Args:
        fusion_power_mw: Total fusion power in megawatt.
        alpha_power_mw: Alpha power deposited in the plasma in megawatt.
        neutron_power_mw: Neutron power leaving the plasma in megawatt.
        auxiliary_power_mw: Heating and current drive power coupled to the
            plasma in megawatt.
        parameters: Engineering parameters.

    Returns:
        The full electrical accounting.

    Raises:
        ValueError: If any input power is negative.
    """
    for name, value in (
        ("fusion_power_mw", fusion_power_mw),
        ("alpha_power_mw", alpha_power_mw),
        ("neutron_power_mw", neutron_power_mw),
        ("auxiliary_power_mw", auxiliary_power_mw),
    ):
        if value < 0.0:
            raise ValueError(f"{name} must be non-negative, got {value}")

    thermal = (
        parameters.blanket_multiplication * neutron_power_mw
        + alpha_power_mw
        + auxiliary_power_mw
    )
    gross = parameters.thermal_efficiency * thermal

    heating_draw = auxiliary_power_mw / parameters.heating_wallplug_efficiency
    balance_of_plant = parameters.balance_of_plant_fraction * gross
    recirculating = (
        heating_draw
        + parameters.cryoplant_mw
        + balance_of_plant
        + parameters.tritium_and_auxiliary_mw
    )
    net = gross - recirculating

    return PlantResult(
        thermal_power_mw=thermal,
        gross_electric_mw=gross,
        heating_draw_mw=heating_draw,
        cryoplant_mw=parameters.cryoplant_mw,
        balance_of_plant_mw=balance_of_plant,
        tritium_and_auxiliary_mw=parameters.tritium_and_auxiliary_mw,
        recirculating_mw=recirculating,
        net_electric_mw=net,
        recirculating_fraction=recirculating / gross if gross > 0.0 else float("inf"),
        net_efficiency=net / fusion_power_mw if fusion_power_mw > 0.0 else float("-inf"),
        engineering_gain=gross / recirculating if recirculating > 0.0 else float("inf"),
    )
