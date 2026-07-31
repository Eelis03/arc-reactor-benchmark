"""Structural interfaces the algorithm layer solves against.

The confinement scaling is the one place where genuinely interchangeable
alternatives exist, and it is the choice that dominates the uncertainty of any
projection, so it is the thing that has to be swappable without touching the
solver.

The Protocol is declared here rather than in the model layer on purpose. The
concrete scalings in :mod:`arc_benchmark.model.confinement` satisfy it
structurally and import nothing from this package, which keeps the dependency
arrow pointing from the solver to the physics and never the other way. A caller
can supply a scaling of its own, from a fit this package has never heard of, and
the solver will accept it as long as it exposes these three members.
"""

from __future__ import annotations

from typing import Protocol, runtime_checkable

from arc_benchmark.model.confinement import ConfinementInputs

__all__ = ["ConfinementScaling"]


@runtime_checkable
class ConfinementScaling(Protocol):
    """An energy confinement time scaling that degrades as a power of loss power.

    The power dependence is exposed separately from the confinement time because
    the steady-state balance is a fixed point in the loss power, and a scaling of
    the form ``tau_E = C P**(-alpha)`` lets that fixed point be solved in closed
    form instead of iterated. Requiring the exponent as part of the interface is
    therefore a deliberate restriction: it excludes scalings that are not power
    laws in the loss power, and it buys an operating point that is reproducible
    to the last bit on any machine.
    """

    @property
    def name(self) -> str:
        """Short identifier, carried through into every trace and report."""
        ...

    @property
    def power_degradation(self) -> float:
        """Positive ``alpha`` such that ``tau_E`` is proportional to ``P**(-alpha)``.

        Must lie strictly between zero and one. At one or above the steady-state
        balance has no solution, because increasing the heating power would not
        increase the loss power.
        """
        ...

    def tau_e_at_unit_power(self, inputs: ConfinementInputs) -> float:
        """Confinement time in second at a loss power of one megawatt."""
        ...

    def tau_e(self, inputs: ConfinementInputs) -> float:
        """Confinement time in second at the loss power carried by ``inputs``."""
        ...
