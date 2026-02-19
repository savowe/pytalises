"""Private propagation engine abstractions.

These classes are internal implementation details and may change without notice.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any


class Engine(ABC):
    """Private engine contract used by propagation kernels."""

    name: str
    xp: Any

    def __init__(self, backend: Any):
        self.backend = backend

    def zeros(self, shape: tuple[int, ...], dtype: Any = "complex128") -> Any:
        return self.backend.zeros(shape, dtype=dtype)

    def asarray(self, data: Any, dtype: Any = None) -> Any:
        return self.backend.asarray(data, dtype=dtype)

    def set_num_threads(self, n: int) -> None:
        self.backend.set_num_threads(n)

    def create_fft_plan(
        self,
        array: Any,
        axes: tuple[int, ...],
        num_threads: int = 1,
        flags: tuple[str, ...] = ("FFTW_ESTIMATE", "FFTW_DESTROY_INPUT"),
    ) -> Any:
        return self.backend.create_fft_plan(
            array,
            axes=axes,
            num_threads=num_threads,
            flags=flags,
        )

    def apply_kinetic_phase(
        self,
        amp: Any,
        *,
        kmesh: tuple[Any, Any, Any],
        alpha: float,
        dt: float,
    ) -> None:
        """Apply kinetic split-step phase in reciprocal space in-place."""
        kx, ky, kz = kmesh
        imag_unit = self.xp.asarray(1j, dtype=amp.dtype)
        phase = self.xp.exp(-imag_unit * alpha * dt * (kx**2 + ky**2 + kz**2))
        amp *= phase[..., self.xp.newaxis]

    def apply_diagonal_phase(self, amp: Any, *, diagonal: Any, dt: float) -> None:
        """Apply diagonal potential phase in-place."""
        imag_unit = self.xp.asarray(1j, dtype=amp.dtype)
        amp *= self.xp.exp(-imag_unit * diagonal * dt)

    def apply_coupled_phase(
        self,
        amp: Any,
        *,
        eigvals: Any,
        eigvecs: Any,
        dt: float,
    ) -> None:
        """Apply non-diagonal potential step via eigendecomposition in-place."""
        # psi' = U @ diag(exp(-i*lambda*dt)) @ U† @ psi
        u_dag = self.xp.swapaxes(self.xp.conjugate(eigvecs), -1, -2)
        tmp = self.xp.matmul(u_dag, amp[..., self.xp.newaxis])[..., 0]
        imag_unit = self.xp.asarray(1j, dtype=amp.dtype)
        tmp *= self.xp.exp(-imag_unit * eigvals * dt)
        amp[...] = self.xp.matmul(eigvecs, tmp[..., self.xp.newaxis])[..., 0]

    def apply_coupled_phase_2x2(
        self,
        amp: Any,
        *,
        matrix: Any,
        dt: float,
        kernel: str = "vectorized",
    ) -> None:
        """Apply non-diagonal potential step for 2x2 Hermitian matrices.

        Uses the closed-form matrix exponential of a Hermitian 2x2 block:
        ``V = c I + [[z, b], [conj(b), -z]]``.
        """
        if amp.shape[-1] != 2 or matrix.shape[-2:] != (2, 2):
            raise ValueError("apply_coupled_phase_2x2 requires 2x2 coupled amplitudes")
        if kernel not in {"vectorized", "fused"}:
            raise ValueError(
                "apply_coupled_phase_2x2 kernel must be 'vectorized' or 'fused'"
            )
        # Base implementation is the vectorized path used by NumPy and as fallback.

        xp = self.xp

        a = matrix[..., 0, 0]
        d = matrix[..., 1, 1]
        b = matrix[..., 0, 1]

        c = 0.5 * (a + d)
        z = 0.5 * (a - d)
        r = xp.sqrt((z * xp.conjugate(z)).real + (b * xp.conjugate(b)).real)

        imag_unit = xp.asarray(1j, dtype=amp.dtype)
        dt_scalar = xp.asarray(dt, dtype=r.dtype)

        phase = xp.exp(-imag_unit * c * dt_scalar)
        cos_term = xp.cos(r * dt_scalar)

        eps = xp.asarray(1e-15, dtype=r.dtype)
        safe_r = xp.where(r > eps, r, xp.asarray(1.0, dtype=r.dtype))
        sin_over_r = xp.sin(r * dt_scalar) / safe_r
        sin_over_r = xp.where(r > eps, sin_over_r, dt_scalar)

        u00 = cos_term - imag_unit * sin_over_r * z
        u11 = cos_term + imag_unit * sin_over_r * z
        u01 = -imag_unit * sin_over_r * b
        u10 = -imag_unit * sin_over_r * xp.conjugate(b)

        psi0 = amp[..., 0]
        psi1 = amp[..., 1]
        new0 = phase * (u00 * psi0 + u01 * psi1)
        new1 = phase * (u10 * psi0 + u11 * psi1)

        amp[..., 0] = new0
        amp[..., 1] = new1

    def inner_product(self, lhs: Any, rhs: Any, *, volume_element: float) -> Any:
        """Return <lhs|rhs> including spatial volume element."""
        return self.xp.sum(lhs * self.xp.conjugate(rhs)) * volume_element

    def state_occupation(self, amp: Any, *, state_index: int, volume_element: float) -> Any:
        """Return occupation of a single internal state."""
        return (
            self.xp.sum(self.xp.abs(amp[:, :, :, state_index]) ** 2)
            * volume_element
        )

    def synchronize(self) -> None:
        """Synchronize outstanding backend work when needed for profiling."""
        return None

    @abstractmethod
    def eigendecompose_hermitian(self, matrices: Any) -> tuple[Any, Any]:
        """Return eigenvalues and eigenvectors for Hermitian blocks."""
