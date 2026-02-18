"""GPU propagation engine."""

from __future__ import annotations

try:
    import cupy as cp
except Exception as exc:  # pragma: no cover - optional dependency
    raise ImportError("GpuEngine requires CuPy.") from exc

from .base import Engine


class GpuEngine(Engine):
    """GPU engine backed by CuPy kernels and linear algebra."""

    name = "gpu"
    xp = cp

    def __init__(self, backend):
        super().__init__(backend)
        device_count = cp.cuda.runtime.getDeviceCount()
        if device_count < 1:
            raise RuntimeError("CuPy backend selected but no CUDA device is available.")

    def eigendecompose_hermitian(self, matrices):
        eigvals, eigvecs = cp.linalg.eigh(matrices)
        return eigvals, eigvecs
