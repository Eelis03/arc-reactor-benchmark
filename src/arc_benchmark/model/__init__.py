"""Physical model: constants, geometry, reactivity, radiation, scalings, limits.

Everything in this subpackage is a pure function or a frozen dataclass. Nothing
here reads a file, writes a figure, or solves a balance. The dependency direction
inside the subpackage is one way: ``constants`` and ``geometry`` know nothing of
anything else, ``reactivity`` uses ``constants``, ``profiles`` uses
``reactivity``, and ``limits`` uses ``constants`` and ``geometry``.
"""

from __future__ import annotations

from arc_benchmark.model.confinement import (
    CONFINEMENT_SCALINGS,
    IPB98Y2,
    ITER89P,
    PETTY08,
    ConfinementInputs,
    PowerLawConfinement,
)
from arc_benchmark.model.geometry import PlasmaGeometry
from arc_benchmark.model.plant import PlantParameters, PlantResult, evaluate_plant
from arc_benchmark.model.profiles import FLAT_PROFILE, ProfileShape
from arc_benchmark.model.radiation import ImpurityRadiator
from arc_benchmark.model.reactivity import cross_section_barn, dt_reactivity

__all__ = [
    "CONFINEMENT_SCALINGS",
    "FLAT_PROFILE",
    "IPB98Y2",
    "ITER89P",
    "PETTY08",
    "ConfinementInputs",
    "ImpurityRadiator",
    "PlantParameters",
    "PlantResult",
    "PlasmaGeometry",
    "PowerLawConfinement",
    "ProfileShape",
    "cross_section_barn",
    "dt_reactivity",
    "evaluate_plant",
]
