"""CPU propagation engine."""

from __future__ import annotations

import numpy as np

from .base import Engine


class CpuEngine(Engine):
    """CPU engine backed by NumPy/pyFFTW backend primitives."""

    name = "cpu"
    xp = np

    def eigendecompose_hermitian(self, matrices):
        return self.backend.eigh(matrices)
