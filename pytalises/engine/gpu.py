"""GPU propagation engine."""

from __future__ import annotations

try:
    import cupy as cp
except Exception as exc:  # pragma: no cover - optional dependency
    raise ImportError("GpuEngine requires CuPy.") from exc

from .base import Engine


@cp.fuse()
def _analytic_2x2_step_fused(psi0, psi1, a, d, b, dt, eps, imag_unit):
    c = 0.5 * (a + d)
    z = 0.5 * (a - d)
    r = cp.sqrt(cp.real(z * cp.conj(z) + b * cp.conj(b)))

    phase = cp.exp(-imag_unit * c * dt)
    cos_term = cp.cos(r * dt)

    safe_r = cp.where(r > eps, r, (eps * 0) + 1)
    sin_over_r = cp.sin(r * dt) / safe_r
    sin_over_r = cp.where(r > eps, sin_over_r, dt)

    u00 = cos_term - imag_unit * sin_over_r * z
    u11 = cos_term + imag_unit * sin_over_r * z
    u01 = -imag_unit * sin_over_r * b
    u10 = -imag_unit * sin_over_r * cp.conj(b)

    new0 = phase * (u00 * psi0 + u01 * psi1)
    new1 = phase * (u10 * psi0 + u11 * psi1)
    return new0, new1


class GpuEngine(Engine):
    """GPU engine backed by CuPy kernels and linear algebra."""

    name = "gpu"
    xp = cp

    def __init__(self, backend):
        super().__init__(backend)
        device_count = cp.cuda.runtime.getDeviceCount()
        if device_count < 1:
            raise RuntimeError("CuPy backend selected but no CUDA device is available.")

    def synchronize(self) -> None:
        cp.cuda.Stream.null.synchronize()

    def eigendecompose_hermitian(self, matrices):
        eigvals, eigvecs = cp.linalg.eigh(matrices)
        return eigvals, eigvecs

    def apply_coupled_phase_2x2(
        self,
        amp,
        *,
        matrix,
        dt: float,
        kernel: str = "vectorized",
    ) -> None:
        if kernel != "fused":
            super().apply_coupled_phase_2x2(amp, matrix=matrix, dt=dt, kernel=kernel)
            return

        if amp.shape[-1] != 2 or matrix.shape[-2:] != (2, 2):
            raise ValueError("apply_coupled_phase_2x2 requires 2x2 coupled amplitudes")

        a = matrix[..., 0, 0]
        d = matrix[..., 1, 1]
        b = matrix[..., 0, 1]
        psi0 = amp[..., 0]
        psi1 = amp[..., 1]

        dt_scalar = cp.asarray(dt, dtype=a.real.dtype)
        eps = cp.asarray(1e-15, dtype=a.real.dtype)
        imag_unit = cp.asarray(1j, dtype=amp.dtype)
        new0, new1 = _analytic_2x2_step_fused(
            psi0, psi1, a, d, b, dt_scalar, eps, imag_unit
        )

        amp[..., 0] = new0
        amp[..., 1] = new1
