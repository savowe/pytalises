"""CuPy backend implementation for pyTALISES.

Notes
-----
This backend requires CuPy and a working CUDA runtime.
"""

from __future__ import annotations

from typing import Any

try:
    import cupy as cp
except ImportError as exc:  # pragma: no cover - exercised via conditional import
    raise ImportError("CupyBackend requires CuPy to be installed.") from exc

from .base import Backend, FFTPlan


class CupyFFTPlan(FFTPlan):
    """Simple FFT plan wrapper using CuPy FFT routines."""

    def __init__(self, array: Any, axes: tuple[int, ...]):
        self._array = array
        self._axes = axes

    def forward(self) -> None:
        self._array[...] = cp.fft.fftn(self._array, axes=self._axes)

    def backward(self) -> None:
        self._array[...] = cp.fft.ifftn(self._array, axes=self._axes)


class CupyBackend(Backend):
    """CuPy backend for GPU-resident propagation."""

    name = "cupy"

    def __init__(self, num_threads: int = 1, device_id: int | None = None):
        self.num_threads = num_threads
        self.device_id = device_id
        if device_id is not None:
            cp.cuda.Device(device_id).use()
        self._eval_globals = {
            "cp": cp,
            "np": cp,
            "pi": cp.pi,
            "exp": cp.exp,
            "sin": cp.sin,
            "cos": cp.cos,
            "tan": cp.tan,
            "arcsin": cp.arcsin,
            "arccos": cp.arccos,
            "arctan": cp.arctan,
            "sinh": cp.sinh,
            "cosh": cp.cosh,
            "tanh": cp.tanh,
            "sqrt": cp.sqrt,
            "log": cp.log,
            "log10": cp.log10,
            "abs": cp.abs,
            "real": cp.real,
            "imag": cp.imag,
            "conjugate": cp.conjugate,
            "conj": cp.conj,
            "where": cp.where,
            "minimum": cp.minimum,
            "maximum": cp.maximum,
            "power": cp.power,
        }

    def set_num_threads(self, n: int) -> None:
        # CuPy/CUDA thread/block scheduling is handled by CUDA runtime.
        self.num_threads = n

    def zeros(self, shape, dtype: Any = "complex128"):
        return cp.zeros(shape, dtype=dtype)

    def empty_aligned(self, shape, dtype: Any = "complex128"):
        return cp.empty(shape, dtype=dtype)

    def asarray(self, data, dtype: Any = None):
        return cp.asarray(data, dtype=dtype)

    def create_fft_plan(
        self,
        array,
        axes,
        num_threads: int = 1,
        flags: tuple[str, ...] = ("FFTW_ESTIMATE", "FFTW_DESTROY_INPUT"),
    ):
        del num_threads, flags
        return CupyFFTPlan(array, axes)

    def exp(self, x):
        return cp.exp(x)

    def einsum(self, subscripts: str, *operands, out=None):
        return cp.einsum(subscripts, *operands, out=out, optimize="greedy")

    def evaluate(self, expr: str, local_dict: dict, global_dict: dict | None = None):
        scope = dict(self._eval_globals)
        if global_dict:
            scope.update(global_dict)
        scope.update(local_dict)
        return eval(expr, {"__builtins__": {}}, scope)

    def eigh(self, matrices):
        eigvals, eigvecs = cp.linalg.eigh(matrices)
        matrices[...] = eigvecs
        return eigvals, matrices

    def meshgrid(self, *arrays, indexing="ij"):
        return cp.meshgrid(*arrays, indexing=indexing)

    def linspace(self, start, stop, num):
        return cp.linspace(start, stop, num)

    def fftfreq(self, n):
        return cp.fft.fftfreq(n)

    def sum(self, x, axis=None):
        return cp.sum(x, axis=axis)

    def abs(self, x):
        return cp.abs(x)

    def sqrt(self, x):
        return cp.sqrt(x)

    def conjugate(self, x):
        return cp.conjugate(x)

    def power(self, x, n):
        return cp.power(x, n)
