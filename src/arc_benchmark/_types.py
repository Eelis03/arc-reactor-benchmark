"""The array alias used by every public signature in this package."""

from __future__ import annotations

import numpy as np
import numpy.typing as npt

__all__ = ["FloatArray"]

FloatArray = npt.NDArray[np.float64]
