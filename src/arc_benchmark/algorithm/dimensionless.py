"""A confinement scaling rewritten in dimensionless plasma parameters.

An empirical scaling is fitted in engineering variables, but the transport it
stands in for is a statement about dimensionless ones. Rewriting the
Bohm-normalised confinement time as a power law in the normalised gyroradius,
the beta, the collisionality, and the safety factor is what makes an empirical
fit comparable with a transport model: Bohm diffusion gives a normalised
gyroradius exponent of minus two and gyro-Bohm diffusion minus three, and a fit
constrained to be independent of beta has a beta exponent of zero by
construction. This package describes Petty08 as exactly that kind of fit, and
this is what checks the description.

The conversion is not a property of the scaling on its own, which is why it
lives in this layer rather than beside the scalings. The loss power in a scaling
is not an independent variable at a solved steady state: it is the stored energy
divided by the confinement time, which is the balance
:mod:`arc_benchmark.algorithm.operating` inverts in closed form. Eliminating it
uses that balance and nothing else.

The derivation runs in two steps, both at fixed inverse aspect ratio,
elongation, and isotope mass, so that those exponents are constants and drop out
of every ratio. Write the scaling as

    tau_E = C I^a_I B^a_B n^a_n P^a_P R^a_R

and substitute ``P = W / tau_E``, with the stored energy ``W`` going as
``n T R**3``. Solving for the confinement time gives the steady-state form

    tau_E = C' I^p_I B^p_B n^p_n T^p_T R^p_R

with every exponent divided by ``1 + a_P``, the density carrying ``a_n + a_P``,
the temperature ``a_P``, and the major radius ``a_R + 3 a_P``. Then, with

    rho* goes as sqrt(T) / (B R),   beta goes as n T / B**2,
    nu* goes as n R / T**2,         q goes as R B / I

matching the exponents of the five engineering variables in

    Omega_ci tau_E, which goes as B tau_E,
        = K rho*^x_rho beta^x_beta nu*^x_nu q^x_q R^x_R

gives five equations in five unknowns with exactly one solution. The last of
them is the exponent of a dimensional length, which no dimensionless parameter
can absorb. It is zero exactly when the scaling satisfies the Kadomtsev
constraint, and its size is how far the fit is from being expressible as
physics rather than as a regression.

Reference:
    ITER Physics Basis, Chapter 2, Nuclear Fusion 39 (1999) 2175, section 6,
    which quotes IPB98(y,2) in this form, and B. B. Kadomtsev, "Tokamaks and
    dimensional analysis", Soviet Journal of Plasma Physics 1 (1975) 295, for
    the constraint the residual measures.
"""

from __future__ import annotations

from dataclasses import dataclass

from arc_benchmark.model.confinement import PowerLawConfinement

__all__ = ["DimensionlessExponents", "dimensionless_exponents"]


@dataclass(frozen=True, slots=True)
class _SteadyStateExponents:
    """The scaling with the loss power eliminated by the steady-state balance."""

    current: float
    field: float
    density: float
    temperature: float
    major_radius: float


@dataclass(frozen=True, slots=True)
class DimensionlessExponents:
    """A confinement scaling as a power law in dimensionless parameters.

    The exponents are those of ``Omega_ci tau_E`` written as
    ``K rho*^x_rho beta^x_beta nu*^x_nu q^x_q R^x_R``, at fixed inverse aspect
    ratio, elongation, and isotope mass. The coefficient ``K`` is not carried:
    it absorbs every constant of proportionality dropped along the way and says
    nothing the exponents do not.

    Attributes:
        scaling_name: Name of the scaling these exponents describe.
        normalised_gyroradius: Exponent of ``rho*``. Minus two is Bohm
            diffusion and minus three is gyro-Bohm, so this is the number that
            says which of the two an empirical fit resembles.
        beta: Exponent of the plasma beta. Zero for a fit constrained to be
            beta-independent, and strongly negative for one that is not.
        collisionality: Exponent of ``nu*``.
        safety_factor: Exponent of the safety factor.
        dimensional_residual: Exponent of the major radius left over once every
            dimensionless parameter has taken its share. Zero exactly when the
            scaling satisfies the Kadomtsev constraint. It is reported rather
            than dropped because its size is the only measure of how far a fit
            is from being expressible in dimensionless variables at all.
    """

    scaling_name: str
    normalised_gyroradius: float
    beta: float
    collisionality: float
    safety_factor: float
    dimensional_residual: float


def _steady_state_exponents(scaling: PowerLawConfinement) -> _SteadyStateExponents:
    """Eliminate the loss power from a scaling using the steady-state balance.

    The stored energy carries one power of density, one of temperature, and
    three of the major radius at fixed shape, so substituting ``P = W / tau_E``
    moves the power degradation onto each of those and divides every exponent by
    ``1 + a_P``.
    """
    alpha = scaling.power_degradation
    if not 0.0 < alpha < 1.0:
        raise ValueError(
            f"scaling {scaling.name!r} has power degradation {alpha}, which must lie strictly "
            "in (0, 1) for the loss power to be eliminated by the steady-state balance"
        )

    e = scaling.exponents
    denominator = 1.0 + e.power
    return _SteadyStateExponents(
        current=e.current / denominator,
        field=e.field / denominator,
        density=(e.density + e.power) / denominator,
        temperature=e.power / denominator,
        major_radius=(e.major_radius + 3.0 * e.power) / denominator,
    )


def dimensionless_exponents(scaling: PowerLawConfinement) -> DimensionlessExponents:
    """Rewrite a confinement scaling in dimensionless plasma parameters.

    The five matching equations are solved in the order they close. The plasma
    current appears in the safety factor and nowhere else, so its exponent comes
    out first. The field equation then fixes one combination of the gyroradius
    and beta exponents, the density equation fixes another of the beta and
    collisionality exponents, and the temperature equation closes the system.
    Whatever the major radius is left with is the residual.

    Args:
        scaling: The scaling to convert. Only its exponents are used, since a
            change of coefficient moves ``K`` and no exponent.

    Returns:
        The exponents of ``rho*``, beta, ``nu*``, and the safety factor,
        together with the residual exponent of the major radius.

    Raises:
        ValueError: If the power degradation is not strictly inside ``(0, 1)``,
            in which case there is no steady state to eliminate the loss power
            with and the conversion is not defined.
    """
    p = _steady_state_exponents(scaling)

    safety_factor = -p.current
    field_equation = 1.0 + p.field + p.current  # equals -(x_rho + 2 x_beta)
    beta = 0.5 * (p.temperature + 2.0 * p.density + 0.5 * field_equation)
    normalised_gyroradius = -field_equation - 2.0 * beta
    collisionality = p.density - beta
    dimensional_residual = (
        p.major_radius + normalised_gyroradius - collisionality - safety_factor
    )

    return DimensionlessExponents(
        scaling_name=scaling.name,
        normalised_gyroradius=normalised_gyroradius,
        beta=beta,
        collisionality=collisionality,
        safety_factor=safety_factor,
        dimensional_residual=dimensional_residual,
    )
