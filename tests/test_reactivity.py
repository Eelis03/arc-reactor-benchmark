"""Tier one: the deuterium tritium reactivity against external references.

Two independent checks are applied. The first compares against values tabulated
in Bosch and Hale (1992), which is external to this package entirely. The second
compares the closed-form reactivity fit against a numerical Maxwellian average of
the separate cross-section fit from the same paper, which is external to this
module: the two fits were derived from the same evaluated data by different
procedures and have different coefficients, so agreement between them tests the
transcription of both.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from arc_benchmark.model.reactivity import (
    REACTIVITY_TEMPERATURE_RANGE_KEV,
    cross_section_barn,
    dt_reactivity,
    dt_reactivity_from_cross_section,
)

_CM3_PER_M3 = 1.0e6

# Bosch and Hale (1992), Table VIII, T(d,n)4He Maxwellian reactivity in cubic
# centimetre per second. The table is printed to four significant figures.
_TABULATED_CM3_S: dict[float, float] = {
    1.0: 6.857e-21,
    2.0: 2.977e-19,
    10.0: 1.136e-16,
    20.0: 4.330e-16,
    50.0: 8.649e-16,
}


@pytest.mark.parametrize(("temperature_kev", "published"), sorted(_TABULATED_CM3_S.items()))
def test_reactivity_matches_the_published_table(temperature_kev: float, published: float) -> None:
    """Reproduce the tabulated reactivity across four orders of magnitude.

    Tolerance: the table is printed to four significant figures, so a value read
    from it carries a rounding uncertainty of half a unit in the fourth digit,
    which is at most 5e-4 in relative terms. The tolerance is set at 1e-3, twice
    that, which leaves room for one digit of transcription slack and is still
    tight enough that a wrong coefficient in the fit could not pass.
    """
    computed = float(dt_reactivity(temperature_kev)) * _CM3_PER_M3
    assert computed == pytest.approx(published, rel=1.0e-3)


def test_reactivity_peaks_where_the_literature_places_it() -> None:
    """The Maxwellian reactivity peaks near 65 keV at about 8.9e-16 cm^3/s.

    Tolerance: the peak location is quoted in the literature to two significant
    figures, so it is checked to within 5 keV, and the peak value to within
    2 percent.
    """
    temperatures = np.linspace(20.0, 100.0, 8001)
    values = dt_reactivity(temperatures)
    peak = int(np.argmax(values))
    assert temperatures[peak] == pytest.approx(65.0, abs=5.0)
    assert float(values[peak]) * _CM3_PER_M3 == pytest.approx(8.9e-16, rel=2.0e-2)


def test_cross_section_resonance_is_at_the_accepted_place() -> None:
    """The cross-section peaks at about 5 barn near 64 keV centre-of-mass energy.

    This is the check that fixes the unit of the cross-section fit. The published
    astrophysical S-function is tabulated in keV millibarn, so an implementation
    that forgets to convert is out by exactly a thousand, which this assertion
    catches immediately.

    Tolerance: the accepted peak is quoted as 5 barn near 64 keV, both to two
    significant figures, so the value is checked to 5 percent and the location to
    2 keV.
    """
    energies = np.linspace(20.0, 200.0, 20001)
    values = cross_section_barn(energies)
    peak = int(np.argmax(values))
    assert float(values[peak]) == pytest.approx(5.0, rel=5.0e-2)
    assert energies[peak] == pytest.approx(64.0, abs=2.0)


@pytest.mark.parametrize("temperature_kev", [1.0, 2.0, 5.0, 10.0, 20.0, 40.0, 70.0, 100.0])
def test_two_independent_fits_agree(temperature_kev: float) -> None:
    """The reactivity fit agrees with the Maxwellian average of the cross-section fit.

    Tolerance: Bosch and Hale state a maximum relative error of 0.25 percent for
    the reactivity fit and about 2 percent for the deuterium tritium
    cross-section fit against the same evaluation, so two fits carrying those
    errors can differ by roughly 2.3 percent without either being wrong.
    Truncating the energy integral at the top of the published cross-section
    range costs a further 0.3 percent at 100 keV. The tolerance is set at
    2.5 percent, the sum of those. The largest deviation actually observed across
    this parameter set is 0.75 percent, so the test has ample margin and would
    still fail on a mistyped coefficient, which moves the answer by far more.
    """
    from_fit = float(dt_reactivity(temperature_kev))
    from_integral = dt_reactivity_from_cross_section(temperature_kev)
    assert from_integral == pytest.approx(from_fit, rel=2.5e-2)


def test_reactivity_is_monotone_below_the_peak() -> None:
    """The reactivity rises without interruption from 1 keV to 40 keV."""
    temperatures = np.linspace(1.0, 40.0, 400)
    values = dt_reactivity(temperatures)
    assert np.all(np.diff(values) > 0.0)


def test_reactivity_rejects_temperatures_outside_the_published_range() -> None:
    """Evaluating outside the fit range raises rather than extrapolating silently."""
    low, high = REACTIVITY_TEMPERATURE_RANGE_KEV
    with pytest.raises(ValueError, match="Bosch and Hale"):
        dt_reactivity(low / 2.0)
    with pytest.raises(ValueError, match="Bosch and Hale"):
        dt_reactivity(high * 2.0)


def test_cross_section_rejects_energies_outside_the_published_range() -> None:
    """The cross-section fit refuses to extrapolate as well."""
    with pytest.raises(ValueError, match="Bosch and Hale"):
        cross_section_barn(0.1)
    with pytest.raises(ValueError, match="Bosch and Hale"):
        cross_section_barn(1000.0)


def test_reactivity_accepts_arrays_and_scalars_alike() -> None:
    """A scalar and a one-element array give the same number."""
    scalar = float(dt_reactivity(14.0))
    array = dt_reactivity(np.array([14.0]))
    assert array.shape == (1,)
    assert float(array[0]) == scalar
    assert math.isfinite(scalar)
