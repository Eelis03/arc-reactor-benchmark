"""Operational limits, written as pure functions of plasma parameters.

Four limits are covered: the Greenwald density limit, the Troyon beta limit, the
edge safety factor, and the bootstrap current fraction. Each is a published
scaling and each is evaluated here without any judgement about whether it is
satisfied. Turning these numbers into verdicts is the job of
:mod:`arc_benchmark.algorithm.constraints`, which keeps the physics separable
from the policy about what counts as an acceptable margin.
"""

from __future__ import annotations

import math
from typing import Final

from arc_benchmark.model.constants import VACUUM_PERMEABILITY
from arc_benchmark.model.geometry import PlasmaGeometry

__all__ = [
    "DEFAULT_BOOTSTRAP_COEFFICIENT",
    "TROYON_COEFFICIENT",
    "average_poloidal_field",
    "bootstrap_fraction",
    "current_at_fixed_q",
    "cylindrical_safety_factor",
    "greenwald_density",
    "lh_threshold_power",
    "normalised_beta",
    "poloidal_beta",
    "safety_factor_95",
    "toroidal_beta",
    "troyon_beta_limit",
]

TROYON_COEFFICIENT: Final[float] = 2.8
"""The original Troyon normalised beta limit in percent metre tesla per megaampere.

Troyon and colleagues (1984) found ideal magnetohydrodynamic stability to fail
above ``beta[%] = g I_p / (a B_0)`` with ``g`` near 2.8 for the equilibria they
scanned. Values between 2.5 and 4 appear in later work depending on the assumed
current and pressure profiles and on whether a conducting wall is present, so the
limit is a parameter of the constraint evaluation rather than a constant of
nature. ARC quotes a design point at 2.59 and ITER at 1.8.
"""

DEFAULT_BOOTSTRAP_COEFFICIENT: Final[float] = 0.7
"""Coefficient of the ``sqrt(epsilon) beta_p`` bootstrap estimate.

The bootstrap fraction depends on the density and temperature profile shapes
through a coefficient that published estimates place between roughly 0.5 and 1.0
for conventional profiles. 0.7 is the middle of that range. The sensitivity of
any conclusion to this number is large and is reported rather than buried.
"""


def greenwald_density(plasma_current_ma: float, minor_radius: float) -> float:
    """Greenwald density limit in inverse cubic metre.

    ``n_G[1e20 m^-3] = I_p[MA] / (pi a[m]**2)``. The limit is an empirical
    boundary above which tokamak discharges disrupt, and it is remarkable for
    depending on the current density and on nothing else: not on field, not on
    heating power, and not on machine size except through the area.

    Args:
        plasma_current_ma: Plasma current in megaampere.
        minor_radius: Minor radius in metre.

    Returns:
        Density limit in inverse cubic metre.

    Reference:
        M. Greenwald et al., "A new look at density limits in tokamaks",
        Nuclear Fusion 28 (1988) 2199.
    """
    if plasma_current_ma <= 0.0:
        raise ValueError(f"plasma_current_ma must be positive, got {plasma_current_ma}")
    if minor_radius <= 0.0:
        raise ValueError(f"minor_radius must be positive, got {minor_radius}")
    return 1.0e20 * plasma_current_ma / (math.pi * minor_radius**2)


def toroidal_beta(pressure_pa: float, toroidal_field: float) -> float:
    """Ratio of plasma pressure to toroidal magnetic pressure, as a fraction.

    ``beta_t = 2 mu_0 p / B_0**2``. Returned as a fraction, not a percentage;
    the percentage convention appears only where the Troyon limit is written.
    """
    if toroidal_field <= 0.0:
        raise ValueError(f"toroidal_field must be positive, got {toroidal_field}")
    return 2.0 * VACUUM_PERMEABILITY * pressure_pa / toroidal_field**2


def poloidal_beta(pressure_pa: float, poloidal_field: float) -> float:
    """Ratio of plasma pressure to poloidal magnetic pressure, as a fraction."""
    if poloidal_field <= 0.0:
        raise ValueError(f"poloidal_field must be positive, got {poloidal_field}")
    return 2.0 * VACUUM_PERMEABILITY * pressure_pa / poloidal_field**2


def average_poloidal_field(geometry: PlasmaGeometry, plasma_current_ma: float) -> float:
    """Poloidal field averaged over the plasma boundary, in tesla.

    ``B_p = mu_0 I_p / L_p`` with ``L_p`` the poloidal circumference. Using the
    true perimeter of the elongated cross-section rather than ``2 pi a`` matters:
    at an elongation of 1.84 the perimeter is 45 percent larger, and the poloidal
    beta that follows from it, which the bootstrap estimate is proportional to,
    is more than a factor of two different.
    """
    if plasma_current_ma <= 0.0:
        raise ValueError(f"plasma_current_ma must be positive, got {plasma_current_ma}")
    return VACUUM_PERMEABILITY * plasma_current_ma * 1.0e6 / geometry.poloidal_perimeter


def normalised_beta(
    beta_toroidal: float,
    minor_radius: float,
    toroidal_field: float,
    plasma_current_ma: float,
) -> float:
    """Normalised beta ``beta_N`` in percent metre tesla per megaampere.

    ``beta_N = beta_t[%] a B_0 / I_p``. This is the quantity the Troyon limit
    bounds, and it is the form in which published design points quote their
    stability margin.
    """
    if plasma_current_ma <= 0.0:
        raise ValueError(f"plasma_current_ma must be positive, got {plasma_current_ma}")
    return 100.0 * beta_toroidal * minor_radius * toroidal_field / plasma_current_ma


def troyon_beta_limit(
    minor_radius: float,
    toroidal_field: float,
    plasma_current_ma: float,
    beta_n_limit: float = TROYON_COEFFICIENT,
) -> float:
    """Maximum toroidal beta permitted by the Troyon scaling, as a fraction.

    ``beta_max[%] = g I_p / (a B_0)``. Returned as a fraction so that it can be
    compared directly against :func:`toroidal_beta`.

    Reference:
        F. Troyon et al., "MHD limits to plasma confinement", Plasma Physics and
        Controlled Fusion 26 (1984) 209.
    """
    if minor_radius <= 0.0 or toroidal_field <= 0.0:
        raise ValueError("minor_radius and toroidal_field must be positive")
    if beta_n_limit <= 0.0:
        raise ValueError(f"beta_n_limit must be positive, got {beta_n_limit}")
    return 0.01 * beta_n_limit * plasma_current_ma / (minor_radius * toroidal_field)


def cylindrical_safety_factor(
    geometry: PlasmaGeometry,
    toroidal_field: float,
    plasma_current_ma: float,
) -> float:
    """Cylindrical safety factor ``q*``, dimensionless.

    ``q* = (2 pi a**2 B_0 / (mu_0 R_0 I_p)) (1 + kappa**2) / 2``. This is the
    elongation-corrected large-aspect-ratio limit, and it is the quantity that
    enters the kink stability condition ``q* > 2``.
    """
    if plasma_current_ma <= 0.0:
        raise ValueError(f"plasma_current_ma must be positive, got {plasma_current_ma}")
    circular = (
        2.0
        * math.pi
        * geometry.minor_radius**2
        * toroidal_field
        / (VACUUM_PERMEABILITY * geometry.major_radius * plasma_current_ma * 1.0e6)
    )
    return circular * (1.0 + geometry.elongation**2) / 2.0


def current_at_fixed_q(
    geometry: PlasmaGeometry,
    toroidal_field: float,
    cylindrical_q: float,
) -> float:
    """Plasma current in megaampere that gives a chosen cylindrical safety factor.

    The exact inverse of :func:`cylindrical_safety_factor`. A sweep over field or
    size has to say what it holds fixed, and holding the safety factor fixed is
    the physically meaningful choice, because it keeps the equilibrium at the
    same distance from the kink boundary while the field or the size changes.
    """
    if cylindrical_q <= 0.0:
        raise ValueError(f"cylindrical_q must be positive, got {cylindrical_q}")
    numerator = (
        2.0
        * math.pi
        * geometry.minor_radius**2
        * toroidal_field
        * (1.0 + geometry.elongation**2)
        / 2.0
    )
    return numerator / (VACUUM_PERMEABILITY * geometry.major_radius * cylindrical_q * 1.0e6)


def safety_factor_95(
    geometry: PlasmaGeometry,
    toroidal_field: float,
    plasma_current_ma: float,
) -> float:
    """Safety factor at the 95 percent flux surface, dimensionless.

    Uses the shaping-corrected form of the ITER Physics Design Guidelines,

        q_95 = (5 a**2 B_0 / (R_0 I_p)) f_shape f_aspect
        f_shape = (1 + kappa**2 (1 + 2 delta**2 - 1.2 delta**3)) / 2
        f_aspect = (1.17 - 0.65 eps) / (1 - eps**2)**2

    with ``I_p`` in megaampere. Evaluated at the ITER reference geometry this
    returns 3.00 against the published 3.0, which is the check the test suite
    applies to it.

    Reference:
        N. A. Uckan and the ITER Physics Group, "ITER physics design guidelines:
        1989", ITER Documentation Series No. 10, IAEA (1990), reproduced in the
        ITER Physics Basis, Nuclear Fusion 39 (1999) 2137.
    """
    if plasma_current_ma <= 0.0:
        raise ValueError(f"plasma_current_ma must be positive, got {plasma_current_ma}")
    eps = geometry.inverse_aspect_ratio
    kappa = geometry.elongation
    delta = geometry.triangularity

    cylindrical = (
        5.0
        * geometry.minor_radius**2
        * toroidal_field
        / (geometry.major_radius * plasma_current_ma)
    )
    shape = (1.0 + kappa**2 * (1.0 + 2.0 * delta**2 - 1.2 * delta**3)) / 2.0
    aspect = (1.17 - 0.65 * eps) / (1.0 - eps**2) ** 2
    return cylindrical * shape * aspect


def lh_threshold_power(
    electron_density: float,
    toroidal_field: float,
    surface_area: float,
) -> float:
    """Power that must cross the separatrix to reach H-mode, in megawatt.

    ``P_LH = 0.0488 n_e20**0.717 B_t**0.803 S**0.941`` with the line-averaged
    density in units of 1e20 inverse cubic metre, the field in tesla, and the
    plasma surface area in square metre.

    This limit is what makes an H-mode confinement scaling conditional rather
    than free. A design point whose loss power falls below this threshold is not
    entitled to an H-mode confinement time, and evaluating IPB98(y,2) there
    produces a number with no support.

    Reference:
        Y. R. Martin, T. Takizuka, and the ITPA CDBM H-mode Threshold Database
        Working Group, "Power requirement for accessing the H-mode in ITER",
        Journal of Physics: Conference Series 123 (2008) 012033.
    """
    if electron_density <= 0.0 or toroidal_field <= 0.0 or surface_area <= 0.0:
        raise ValueError("density, field, and surface area must all be positive")
    n_e20 = electron_density / 1.0e20
    return float(0.0488 * n_e20**0.717 * toroidal_field**0.803 * surface_area**0.941)


def bootstrap_fraction(
    geometry: PlasmaGeometry,
    beta_poloidal: float,
    coefficient: float = DEFAULT_BOOTSTRAP_COEFFICIENT,
) -> float:
    """Fraction of the plasma current driven by the bootstrap effect.

    ``f_BS = c sqrt(epsilon) beta_p``. The square root of the inverse aspect
    ratio comes from the trapped particle fraction, and the poloidal beta
    measures how much pressure gradient there is to drive the current.

    The result is not clipped to one. A value above one means the model has been
    handed a plasma whose pressure gradient would drive more current than the
    equilibrium carries, which is a design point to reject, not a number to hide.
    The constraint layer reports it as a violation.

    Reference:
        The form and the coefficient range are those of the ITER Physics Basis,
        Chapter 2, Nuclear Fusion 39 (1999) 2175, section 3.2, following
        H. R. Wilson, "Bootstrap current scaling in tokamaks", Nuclear Fusion 32
        (1992) 257.
    """
    if coefficient <= 0.0:
        raise ValueError(f"coefficient must be positive, got {coefficient}")
    return coefficient * math.sqrt(geometry.inverse_aspect_ratio) * beta_poloidal
