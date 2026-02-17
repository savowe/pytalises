"""Backend abstraction interfaces for pyTALISES."""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any, Tuple


class FFTPlan(ABC):
    """Abstract FFT plan interface."""

    @abstractmethod
    def forward(self) -> None:
        """Execute forward FFT in-place."""

    @abstractmethod
    def backward(self) -> None:
        """Execute backward FFT in-place."""


class Backend(ABC):
    """Abstract backend interface."""

    name: str

    @abstractmethod
    def zeros(self, shape: Tuple[int, ...], dtype: Any = "complex128") -> Any:
        """Create zero-filled array."""

    @abstractmethod
    def empty_aligned(self, shape: Tuple[int, ...], dtype: Any = "complex128") -> Any:
        """Create uninitialized aligned array."""

    @abstractmethod
    def asarray(self, data: Any, dtype: Any = None) -> Any:
        """Convert data to backend-native array."""

    @abstractmethod
    def create_fft_plan(
        self,
        array: Any,
        axes: Tuple[int, ...],
        num_threads: int = 1,
        flags: tuple[str, ...] = ("FFTW_ESTIMATE", "FFTW_DESTROY_INPUT"),
    ) -> FFTPlan:
        """Create forward/backward FFT plan for array."""

    @abstractmethod
    def exp(self, x: Any) -> Any:
        """Element-wise exponential."""

    @abstractmethod
    def einsum(self, subscripts: str, *operands: Any, out: Any = None) -> Any:
        """Einstein summation."""

    @abstractmethod
    def evaluate(self, expr: str, local_dict: dict, global_dict: dict | None = None) -> Any:
        """Evaluate expression with backend expression engine."""

    @abstractmethod
    def eigh(self, matrices: Any) -> tuple[Any, Any]:
        """Batched eigendecomposition of Hermitian matrices."""

    @abstractmethod
    def meshgrid(self, *arrays: Any, indexing: str = "ij") -> list[Any]:
        """Create meshgrid arrays."""

    @abstractmethod
    def linspace(self, start: float, stop: float, num: int) -> Any:
        """Evenly spaced samples."""

    @abstractmethod
    def fftfreq(self, n: int) -> Any:
        """FFT sample frequencies."""

    @abstractmethod
    def sum(self, x: Any, axis: Any = None) -> Any:
        """Sum reduction."""

    @abstractmethod
    def abs(self, x: Any) -> Any:
        """Absolute value."""

    @abstractmethod
    def sqrt(self, x: Any) -> Any:
        """Square root."""

    @abstractmethod
    def conjugate(self, x: Any) -> Any:
        """Complex conjugate."""

    @abstractmethod
    def power(self, x: Any, n: Any) -> Any:
        """Element-wise power."""

    @abstractmethod
    def set_num_threads(self, n: int) -> None:
        """Configure backend threading where relevant."""
