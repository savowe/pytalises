"""NumPy backend implementation for pyTALISES."""

from __future__ import annotations

import numexpr as ne
from numba import jit, prange, set_num_threads
import numpy as np
from numpy.linalg import eigh
import pyfftw

from .base import Backend, FFTPlan


class NumpyFFTPlan(FFTPlan):
    """pyFFTW-based FFT plan wrapper."""

    def __init__(
        self,
        array,
        axes,
        num_threads=1,
        flags=("FFTW_ESTIMATE", "FFTW_DESTROY_INPUT"),
    ):
        self._forward = pyfftw.FFTW(
            array,
            array,
            axes=axes,
            direction="FFTW_FORWARD",
            threads=num_threads,
            flags=flags,
        )
        self._backward = pyfftw.FFTW(
            array,
            array,
            axes=axes,
            direction="FFTW_BACKWARD",
            threads=num_threads,
            flags=flags,
        )

    def forward(self) -> None:
        self._forward()

    def backward(self) -> None:
        self._backward()


class NumpyBackend(Backend):
    """NumPy/NumExpr/pyFFTW backend preserving current behavior."""

    name = "numpy"

    def __init__(self, num_threads: int = 1):
        self.num_threads = num_threads
        self.set_num_threads(num_threads)

    def set_num_threads(self, n: int) -> None:
        self.num_threads = n
        set_num_threads(n)
        ne.set_num_threads(n)

    def zeros(self, shape, dtype="complex128"):
        return np.zeros(shape, dtype=dtype, order="C")

    def empty_aligned(self, shape, dtype="complex128"):
        return pyfftw.empty_aligned(shape, dtype=dtype, order="C")

    def asarray(self, data, dtype=None):
        return np.asarray(data, dtype=dtype)

    def create_fft_plan(
        self,
        array,
        axes,
        num_threads=1,
        flags=("FFTW_ESTIMATE", "FFTW_DESTROY_INPUT"),
    ):
        pyfftw.config.NUM_THREADS = num_threads
        return NumpyFFTPlan(array, axes, num_threads, flags)

    def exp(self, x):
        return np.exp(x)

    def einsum(self, subscripts, *operands, out=None):
        return np.einsum(
            subscripts,
            *operands,
            out=out,
            optimize="optimal",
            order="C",
        )

    def evaluate(self, expr, local_dict, global_dict=None):
        return ne.evaluate(
            expr,
            local_dict=local_dict,
            global_dict=global_dict or {},
            order="C",
        )

    def eigh(self, matrices):
        eigvals = np.zeros(matrices.shape[:-1], order="C", dtype="complex128")
        _batched_eigh_numba(matrices, eigvals)
        return eigvals, matrices

    def meshgrid(self, *arrays, indexing="ij"):
        return np.meshgrid(*arrays, indexing=indexing)

    def linspace(self, start, stop, num):
        return np.linspace(start, stop, num)

    def fftfreq(self, n):
        return np.fft.fftfreq(n)

    def sum(self, x, axis=None):
        return np.sum(x, axis=axis)

    def abs(self, x):
        return np.abs(x)

    def sqrt(self, x):
        return np.sqrt(x)

    def conjugate(self, x):
        return np.conjugate(x)

    def power(self, x, n):
        return np.power(x, n)


@jit(nopython=True, parallel=True, nogil=True, fastmath=True)
def _batched_eigh_numba(matrices, eigvals):
    nX, nY, nZ = matrices.shape[:3]
    for i in prange(nX):
        for j in prange(nY):
            for k in prange(nZ):
                eigvals[i, j, k, :], matrices[i, j, k, :, :] = eigh(
                    matrices[i, j, k, :, :]
                )
