"""Solvers: power balance, operating point, constraint verdicts, Lawson condition.

This layer consumes the model layer and produces solved states. It does no
plotting and writes no files. The confinement scaling reaches it through the
:class:`~arc_benchmark.algorithm.protocols.ConfinementScaling` Protocol, so the
solver never names a particular fit.
"""

from __future__ import annotations

from arc_benchmark.algorithm.balance import (
    LossPowerConvention,
    PlasmaComposition,
    PlasmaState,
    PowerTerms,
    power_terms,
)
from arc_benchmark.algorithm.constraints import (
    ConstraintCheck,
    ConstraintLimits,
    ConstraintReport,
    ConstraintSense,
    cylindrical_q,
    evaluate_constraints,
)
from arc_benchmark.algorithm.lawson import (
    LawsonPoint,
    lawson_n_tau,
    lawson_triple_product,
    optimum_lawson_temperature,
)
from arc_benchmark.algorithm.operating import (
    OperatingPoint,
    solve_ignition_temperature,
    solve_operating_point,
)
from arc_benchmark.algorithm.protocols import ConfinementScaling

__all__ = [
    "ConfinementScaling",
    "ConstraintCheck",
    "ConstraintLimits",
    "ConstraintReport",
    "ConstraintSense",
    "LawsonPoint",
    "LossPowerConvention",
    "OperatingPoint",
    "PlasmaComposition",
    "PlasmaState",
    "PowerTerms",
    "cylindrical_q",
    "evaluate_constraints",
    "lawson_n_tau",
    "lawson_triple_product",
    "optimum_lawson_temperature",
    "power_terms",
    "solve_ignition_temperature",
    "solve_operating_point",
]
