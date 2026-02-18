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
        phase = self.xp.exp(-1j * alpha * dt * (kx**2 + ky**2 + kz**2))
        amp *= phase[..., self.xp.newaxis]

    def apply_diagonal_phase(self, amp: Any, *, diagonal: Any, dt: float) -> None:
        """Apply diagonal potential phase in-place."""
        amp *= self.xp.exp(-1j * diagonal * dt)

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
        tmp *= self.xp.exp(-1j * eigvals * dt)
        amp[...] = self.xp.matmul(eigvecs, tmp[..., self.xp.newaxis])[..., 0]

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
