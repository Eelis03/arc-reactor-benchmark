"""Tier one: the confinement scalings, including an external check against ITER.

The check that matters most here is that IPB98(y,2) evaluated at the ITER
reference parameters returns the confinement time ITER publishes. That number was
produced by the people who fitted the scaling, from the same scaling, so
reproducing it tests the transcription of every exponent and the coefficient at
once.
"""

from __future__ import annotations

import dataclasses

import pytest

from arc_benchmark.algorithm.protocols import ConfinementScaling
from arc_benchmark.model.confinement import (
    CONFINEMENT_SCALINGS,
    IPB98Y2,
    ITER89P,
    PETTY08,
    ConfinementInputs,
)

# The ITER inductive baseline, in the units IPB98(y,2) is published in.
_ITER_REFERENCE = ConfinementInputs(
    plasma_current_ma=15.0,
    toroidal_field=5.3,
    line_averaged_density_e19=10.0,
    loss_power_mw=87.0,
    major_radius=6.2,
    inverse_aspect_ratio=2.0 / 6.2,
    areal_elongation=1.7,
    mass_number=2.5,
)

_ITER_PUBLISHED_TAU_S = 3.7


def test_ipb98_reproduces_the_iter_reference_confinement_time() -> None:
    """IPB98(y,2) at the ITER baseline returns the published 3.7 s.

    Tolerance: the loss power of the ITER baseline is quoted between 80 and
    90 MW across the sources that give it, depending on how much core radiation
    is subtracted. Over that interval this scaling returns 3.81 s down to 3.51 s,
    a spread of plus or minus 4 percent about the published value. The tolerance
    is therefore 5 percent, taken from the uncertainty in the input rather than
    from the error that happens to be observed, which is 2.9 percent.
    """
    tau = IPB98Y2.tau_e(_ITER_REFERENCE)
    assert tau == pytest.approx(_ITER_PUBLISHED_TAU_S, rel=5.0e-2)


def test_the_iter_published_time_lies_inside_the_quoted_loss_power_range() -> None:
    """Sweeping the loss power over its published range brackets 3.7 s.

    This is the statement the tolerance above was derived from, asserted rather
    than only written down, so that the derivation cannot quietly stop being
    true.
    """
    at_80 = IPB98Y2.tau_e(_ITER_REFERENCE.with_loss_power(80.0))
    at_90 = IPB98Y2.tau_e(_ITER_REFERENCE.with_loss_power(90.0))
    assert at_90 < _ITER_PUBLISHED_TAU_S < at_80


@pytest.mark.parametrize("scaling", sorted(CONFINEMENT_SCALINGS.values(), key=lambda s: s.name))
@pytest.mark.parametrize(
    ("field_name", "attribute"),
    [
        ("plasma_current_ma", "current"),
        ("toroidal_field", "field"),
        ("line_averaged_density_e19", "density"),
        ("loss_power_mw", "power"),
        ("major_radius", "major_radius"),
        ("areal_elongation", "elongation"),
        ("mass_number", "mass"),
    ],
)
def test_each_exponent_is_recovered_by_perturbing_its_input(
    scaling: object, field_name: str, attribute: str
) -> None:
    """Doubling one input multiplies the confinement time by two to its exponent.

    This is the unit consistency check applied one variable at a time. It is
    stronger than reading the exponents off the dataclass, because it goes
    through the same arithmetic the solver uses.

    Tolerance: the evaluation is a product of eight powers, so the relative error
    is bounded by roughly sixteen units in the last place. 1e-12 is four orders
    of magnitude above that and still far below any exponent error.
    """
    assert isinstance(scaling, ConfinementScaling)
    exponent = getattr(scaling.exponents, attribute)  # type: ignore[attr-defined]
    doubled = dataclasses.replace(
        _ITER_REFERENCE, **{field_name: 2.0 * getattr(_ITER_REFERENCE, field_name)}
    )
    ratio = scaling.tau_e(doubled) / scaling.tau_e(_ITER_REFERENCE)
    assert ratio == pytest.approx(2.0**exponent, rel=1.0e-12)


def test_inverse_aspect_ratio_exponent_is_recovered() -> None:
    """The aspect ratio exponent, checked separately because it cannot be doubled.

    Doubling the ITER inverse aspect ratio would exceed one, so it is scaled by
    a factor of 1.25 instead and the same identity is applied.
    """
    factor = 1.25
    for scaling in CONFINEMENT_SCALINGS.values():
        perturbed = dataclasses.replace(
            _ITER_REFERENCE,
            inverse_aspect_ratio=factor * _ITER_REFERENCE.inverse_aspect_ratio,
        )
        ratio = scaling.tau_e(perturbed) / scaling.tau_e(_ITER_REFERENCE)
        expected = factor**scaling.exponents.inverse_aspect_ratio
        assert ratio == pytest.approx(expected, rel=1.0e-12)


def test_unit_power_factorisation_is_exact() -> None:
    """``tau_e`` equals ``tau_e_at_unit_power`` times the power factor.

    The operating point solver relies on this identity to invert the balance in
    closed form, so it is asserted rather than assumed.
    """
    for scaling in CONFINEMENT_SCALINGS.values():
        expected = _ITER_REFERENCE.loss_power_mw**-scaling.power_degradation
        ratio = scaling.tau_e(_ITER_REFERENCE) / scaling.tau_e_at_unit_power(_ITER_REFERENCE)
        assert ratio == pytest.approx(expected, rel=1.0e-14)


def test_every_scaling_degrades_with_power_but_not_too_fast() -> None:
    """Power degradation lies strictly inside zero and one for every scaling.

    Outside that interval the steady-state balance has no solution, so this is a
    precondition of the solver rather than a property of the physics.
    """
    for scaling in CONFINEMENT_SCALINGS.values():
        assert 0.0 < scaling.power_degradation < 1.0


def test_l_mode_confinement_is_shorter_than_h_mode() -> None:
    """The L-mode fit returns a shorter time than either H-mode fit at ITER."""
    l_mode = ITER89P.tau_e(_ITER_REFERENCE)
    assert l_mode < IPB98Y2.tau_e(_ITER_REFERENCE)
    assert l_mode < PETTY08.tau_e(_ITER_REFERENCE)


def test_petty_degrades_less_with_power_than_ipb98() -> None:
    """Petty08 has the weaker power degradation, which is why it changes the answer.

    A weaker degradation means the loss power at a solved steady state is less
    sensitive to the stored energy, so the same design point comes out with a
    higher gain under Petty08 than under IPB98(y,2).
    """
    assert PETTY08.power_degradation < IPB98Y2.power_degradation


def test_scalings_satisfy_the_solver_protocol() -> None:
    """Every scaling is accepted by the structural interface the solver needs."""
    for scaling in CONFINEMENT_SCALINGS.values():
        assert isinstance(scaling, ConfinementScaling)


def test_confinement_inputs_reject_unphysical_values() -> None:
    """Every input a power law cannot be evaluated at is refused."""
    with pytest.raises(ValueError, match="plasma_current_ma"):
        dataclasses.replace(_ITER_REFERENCE, plasma_current_ma=0.0)
    with pytest.raises(ValueError, match="loss_power_mw"):
        dataclasses.replace(_ITER_REFERENCE, loss_power_mw=-1.0)
    with pytest.raises(ValueError, match="inverse_aspect_ratio"):
        dataclasses.replace(_ITER_REFERENCE, inverse_aspect_ratio=1.5)
