"""Power balance and confinement benchmarking for ARC-class compact high-field tokamaks.

ARC here is the Affordable, Robust, Compact tokamak concept of Sorbom and
colleagues (2015), a peer-reviewed high-field compact fusion pilot plant design.
It is not the fictional device of a similar nickname.

The package is arranged in five layers with a one-way dependency:
``model`` holds pure physics, ``algorithm`` solves balances against it,
``pipeline`` runs sweeps and benchmark cases, ``analysis`` reads the resulting
traces, and the scripts in ``examples/`` only wire those together.
"""

from __future__ import annotations

from arc_benchmark.algorithm.balance import (
    LossPowerConvention,
    PlasmaComposition,
    PlasmaState,
)
from arc_benchmark.algorithm.constraints import ConstraintLimits, evaluate_constraints
from arc_benchmark.algorithm.operating import OperatingPoint, solve_operating_point
from arc_benchmark.model.confinement import IPB98Y2, ITER89P, PETTY08
from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.model.plant import PlantParameters
from arc_benchmark.model.profiles import ProfileShape
from arc_benchmark.pipeline.machines import machine

__all__ = [
    "IPB98Y2",
    "ITER89P",
    "PETTY08",
    "ConstraintLimits",
    "LossPowerConvention",
    "OperatingPoint",
    "PlantParameters",
    "PlasmaComposition",
    "PlasmaGeometry",
    "PlasmaState",
    "ProfileShape",
    "__version__",
    "evaluate_constraints",
    "machine",
    "solve_operating_point",
]

__version__ = "0.1.0"
