"""v1 compatibility layer for pyTALISES.

This module preserves the pre-v2 call signatures for one migration cycle.
"""

from __future__ import annotations

from typing import Any
import warnings

import numpy as np

from pytalises.grid import Grid
from pytalises.options import PropagationOptions
from pytalises.potentials import (
    BasePotential,
    DiagonalPotential,
    HermitianPotential,
)
from pytalises.propagator import Propagator as _Propagator
from pytalises.propagator import freely_propagate as _freely_propagate
from pytalises.propagator import propagate as _propagate
from pytalises.wavefunction import Wavefunction as _Wavefunction


_LEGACY_MSG = (
    "pytalises.legacy provides compatibility with the pre-v2 API and will be "
    "removed in a future release. Please migrate to the v2 API "
    "(Grid, structured potentials, steps/dt, PropagationOptions)."
)


def _warn_legacy() -> None:
    warnings.warn(_LEGACY_MSG, DeprecationWarning, stacklevel=3)


def _normalize_grid(
    number_of_grid_points: tuple[int, ...] | int,
    spatial_ext: list[tuple[float, float]] | tuple[float, float],
) -> Grid:
    if isinstance(number_of_grid_points, int):
        number_of_grid_points = (number_of_grid_points,)
    shape_input = tuple(number_of_grid_points)

    if isinstance(spatial_ext, tuple) and len(spatial_ext) == 2 and all(
        isinstance(x, (float, int)) for x in spatial_ext
    ):
        ext_input = [tuple(float(v) for v in spatial_ext)]
    else:
        ext_input = [tuple(float(v) for v in pair) for pair in spatial_ext]  # type: ignore[arg-type]

    active_indices = [i for i, n in enumerate(shape_input) if n > 1]
    if not active_indices:
        active_indices = [0]
        if not shape_input:
            shape_input = (1,)

    active_shape = tuple(shape_input[i] for i in active_indices)

    if len(ext_input) == len(shape_input):
        active_extent = tuple(ext_input[i] for i in active_indices)
    elif len(ext_input) == len(active_shape):
        active_extent = tuple(ext_input)
    else:
        raise ValueError(
            "spatial_ext must match either full dimensionality of number_of_grid_points "
            "or only the active (n>1) dimensions."
        )

    return Grid(shape=active_shape, extent=active_extent)


def _to_structured_potential(
    potential: BasePotential | str | list[str],
    *,
    diag: bool,
    num_states: int,
) -> BasePotential:
    if isinstance(potential, BasePotential):
        return potential

    if isinstance(potential, str):
        exprs = [potential]
    else:
        exprs = list(potential)

    if diag:
        return DiagonalPotential(exprs)
    return HermitianPotential.from_lower_triangular(exprs, num_states=num_states)


class Wavefunction(_Wavefunction):
    """Legacy Wavefunction signature wrapper."""

    def __init__(
        self,
        psi: str | list[str],
        number_of_grid_points: tuple[int, ...] | int,
        spatial_ext: list[tuple[float, float]] | tuple[float, float],
        t0: float = 0.0,
        m: float = 1.054571817e-34,
        variables: dict[str, Any] | None = None,
        normalize_const: float | None = None,
    ) -> None:
        _warn_legacy()
        grid = _normalize_grid(number_of_grid_points, spatial_ext)
        super().__init__(
            initial=psi,
            grid=grid,
            t0=t0,
            m=m,
            variables=variables,
            normalize_const=normalize_const,
        )

    @property
    def r(self):  # type: ignore[override]
        """Legacy behavior: array for 1D, list with padded axes otherwise."""
        if self.num_ext_dim == 1:
            return super().r[0]
        return self._r

    @property
    def k(self):  # type: ignore[override]
        """Legacy behavior: array for 1D, list with padded axes otherwise."""
        if self.num_ext_dim == 1:
            return super().k[0]
        return self._k

    @property
    def amp(self):  # type: ignore[override]
        """Legacy behavior: squeezed amplitude view."""
        return np.squeeze(super().amp)

    def freely_propagate(
        self,
        num_time_steps: int,
        delta_t: float,
        num_of_threads: int = 1,
        FFTWflags: tuple[str, ...] = (
            "FFTW_ESTIMATE",
            "FFTW_DESTROY_INPUT",
        ),
    ) -> None:
        options = PropagationOptions(threads=num_of_threads, fftw_flags=FFTWflags)
        super().freely_propagate(steps=num_time_steps, dt=delta_t, options=options)

    def propagate(
        self,
        potential: BasePotential | str | list[str],
        num_time_steps: int,
        delta_t: float,
        **kwargs: Any,
    ) -> None:
        variables = kwargs.get("variables")
        diag = kwargs.get("diag", False)
        num_of_threads = kwargs.get("num_of_threads", 1)
        FFTWflags = kwargs.get(
            "FFTWflags",
            (
                "FFTW_ESTIMATE",
                "FFTW_DESTROY_INPUT",
            ),
        )

        options = PropagationOptions(threads=num_of_threads, fftw_flags=FFTWflags)
        structured = _to_structured_potential(
            potential,
            diag=diag,
            num_states=self.num_int_dim,
        )
        super().propagate(
            potential=structured,
            steps=num_time_steps,
            dt=delta_t,
            variables=variables,
            options=options,
        )


class Propagator:
    """Legacy Propagator signature wrapper."""

    Potential = _Propagator.Potential

    def __init__(
        self,
        psi: _Wavefunction,
        potential: BasePotential | str | list[str],
        variables: dict[str, Any] | None = None,
        diag: bool = False,
        num_of_threads: int = 1,
        FFTWflags: tuple[str, ...] = (
            "FFTW_ESTIMATE",
            "FFTW_DESTROY_INPUT",
        ),
    ) -> None:
        _warn_legacy()
        options = PropagationOptions(threads=num_of_threads, fftw_flags=FFTWflags)
        structured = _to_structured_potential(
            potential,
            diag=diag,
            num_states=psi.num_int_dim,
        )
        self._inner = _Propagator(
            psi,
            structured,
            variables=variables,
            options=options,
        )
        self.psi = self._inner.psi
        self.v = self._inner.v

    def potential_prop(self, delta_t: float) -> None:
        self._inner.potential_prop(delta_t)

    def kinetic_prop(self, delta_t: float) -> None:
        self._inner.kinetic_prop(delta_t)

    def __getattr__(self, name: str):
        return getattr(self._inner, name)


def propagate(
    psi: _Wavefunction,
    potential: BasePotential | str | list[str],
    num_time_steps: int,
    delta_t: float,
    **kwargs: Any,
) -> None:
    """Legacy module-level propagate wrapper."""
    _warn_legacy()
    variables = kwargs.get("variables")
    diag = kwargs.get("diag", False)
    num_of_threads = kwargs.get("num_of_threads", 1)
    FFTWflags = kwargs.get(
        "FFTWflags",
        (
            "FFTW_ESTIMATE",
            "FFTW_DESTROY_INPUT",
        ),
    )
    options = PropagationOptions(threads=num_of_threads, fftw_flags=FFTWflags)
    structured = _to_structured_potential(
        potential,
        diag=diag,
        num_states=psi.num_int_dim,
    )
    _propagate(
        psi,
        structured,
        num_time_steps,
        delta_t,
        variables=variables,
        options=options,
    )


def freely_propagate(
    psi: _Wavefunction,
    num_time_steps: int,
    delta_t: float,
    num_of_threads: int = 1,
    FFTWflags: tuple[str, ...] = (
        "FFTW_ESTIMATE",
        "FFTW_DESTROY_INPUT",
    ),
) -> None:
    """Legacy module-level free propagation wrapper."""
    _warn_legacy()
    options = PropagationOptions(threads=num_of_threads, fftw_flags=FFTWflags)
    _freely_propagate(
        psi,
        num_time_steps,
        delta_t,
        options=options,
    )
