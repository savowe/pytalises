"""The Wavefunction class and its attributes."""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

import pytalises.propagator
from pytalises.backends import get_backend
from pytalises.backends.base import Backend
from pytalises.engine import ExpressionEvaluator, create_engine
from pytalises.grid import Grid
from pytalises.options import PropagationOptions
from pytalises.potentials import BasePotential


class Wavefunction:
    """Wavefunction state for split-step propagation.

    Parameters
    ----------
    initial
        Initial amplitude expression(s). One string per internal state.
    grid
        Explicit spatial grid definition.
    t0
        Initial time.
    m
        Particle mass.
    variables
        Additional variables used in expressions.
    normalize_const
        If provided, normalize total population to this value.
    backend
        Backend name (``"auto"``, ``"numpy"``) or backend instance.
    num_threads
        Thread count for backend routines.
    """

    def __init__(
        self,
        initial: str | list[str],
        grid: Grid,
        t0: float = 0.0,
        m: float = 1.054571817e-34,
        variables: dict[str, Any] | None = None,
        normalize_const: float | None = None,
        backend: str | Backend | None = "auto",
        num_threads: int = 1,
    ) -> None:
        if variables is None:
            variables = {}
        if isinstance(initial, str):
            initial = [initial]
        if not isinstance(grid, Grid):
            raise TypeError("grid must be an instance of pytalises.Grid")
        if not isinstance(t0, float):
            raise TypeError("t0 must be float")
        if not isinstance(m, float):
            raise TypeError("m must be float")
        if not isinstance(variables, dict):
            raise TypeError("variables must be dict")

        self._backend = get_backend(backend, num_threads=num_threads)
        self._engine = create_engine(self._backend)
        self._evaluator = ExpressionEvaluator(backend_name=self._backend.name)

        self.grid = grid
        self.initial = list(initial)
        self.num_int_dim = len(self.initial)

        self.number_of_grid_points = self.grid.padded_shape
        self.spatial_ext = list(self.grid.padded_extent)

        self.axes = tuple(1 if n > 1 else 0 for n in self.number_of_grid_points)
        self.num_ext_dim = int(sum(self.axes))
        self.nX, self.nY, self.nZ = self.number_of_grid_points

        r: list[Any] = []
        Delta_r: list[float] = []
        k: list[Any] = []
        delta_k: list[float] = []

        for i, spatial_ext_tuple in enumerate(self.spatial_ext):
            r_min = spatial_ext_tuple[0]
            r_max = spatial_ext_tuple[1]
            n = self.number_of_grid_points[i]
            r_axis = self._backend.linspace(r_min, r_max, num=n)
            r.append(r_axis)
            if self.axes[i] == 0:
                Delta_r.append(np.nan)
                delta_k.append(np.nan)
                # Keep inactive reciprocal axes as 1D arrays for backend parity
                # (CuPy meshgrid requires array-like inputs, not Python scalars).
                k.append(self._backend.asarray([0.0]))
            else:
                drange = r_max - r_min
                Delta_r.append(drange)
                delta_k.append(2 * np.pi / drange)
                k.append(self._backend.fftfreq(n) * 2 * np.pi * n / drange)

        self._r = r
        self.Delta_r = Delta_r
        self.delta_r = [
            Delta / self.number_of_grid_points[i]
            for i, Delta in enumerate(self.Delta_r)
        ]
        self.delta_k = delta_k
        self._k = k

        self.rmesh = self._backend.meshgrid(*r, indexing="ij")
        self.kmesh = self._backend.meshgrid(*k, indexing="ij")

        self._amp = self._backend.empty_aligned(
            self.number_of_grid_points + (self.num_int_dim,),
            dtype="complex128",
        )

        self.t = t0
        self.m = m
        self.alpha = 1.054571817e-34 / (2 * self.m)

        self.default_var_dict = {
            "alpha": self.alpha,
            "x": self.rmesh[0],
            "y": self.rmesh[1],
            "z": self.rmesh[2],
        }
        for i in range(self.num_int_dim):
            self.default_var_dict[f"psi{i}"] = self._amp[:, :, :, i]

        self.variables = variables
        for i in range(self.num_int_dim):
            self._amp[:, :, :, i] = self._evaluator.eval(
                self.initial[i],
                local_dict={**self.default_var_dict, **self.variables},
                global_dict={"t": self.t},
            )

        self.normalize_const = normalize_const
        if normalize_const is not None:
            self.normalize_to(normalize_const)

        self.construct_FFT(num_threads=num_threads)

    def construct_FFT(
        self,
        num_threads: int = 1,
        FFTWflags: tuple[str, ...] = (
            "FFTW_ESTIMATE",
            "FFTW_DESTROY_INPUT",
        ),
    ) -> None:
        """Construct FFT bindings through engine plan interface."""
        axes = tuple(i for i, active in enumerate(self.axes) if active)
        self._engine.set_num_threads(num_threads)
        self._fft_plan = self._engine.create_fft_plan(
            self._amp,
            axes=axes,
            num_threads=num_threads,
            flags=FFTWflags,
        )
        self.fft = self._fft_plan.forward
        self.ifft = self._fft_plan.backward

    def _active_axis_indices(self) -> list[int]:
        return [i for i, active in enumerate(self.axes) if active]

    def _volume_element(self) -> float:
        elements = [self.delta_r[i] for i, active in enumerate(self.axes) if active]
        if not elements:
            return 1.0
        return float(np.prod(elements))

    @property
    def r(self) -> tuple[NDArray[np.floating[Any]], ...]:
        """Spatial coordinate arrays per active axis."""
        return tuple(self._r[i] for i in self._active_axis_indices())

    @property
    def k(self) -> tuple[NDArray[np.floating[Any]], ...]:
        """Reciprocal coordinate arrays per active axis."""
        return tuple(self._k[i] for i in self._active_axis_indices())

    @property
    def amp(self) -> NDArray[np.complexfloating[Any, Any]]:
        """Wavefunction amplitudes in canonical public shape.

        Shape is ``grid.shape + (num_internal_states,)``.
        """
        ext_dims = self.grid.ndim
        index = (
            (slice(None),) * ext_dims + (0,) * (3 - ext_dims) + (slice(None),)
        )
        return self._amp[index]

    def exp_pos(self, axis: int | None = None) -> NDArray[np.floating[Any]]:
        """Calculate expected position for one axis or all external axes."""
        active_axes = self._active_axis_indices()
        if axis is None:
            exp_pos = np.empty((self.num_ext_dim,))
            for i in range(self.num_ext_dim):
                exp_pos[i] = self.exp_pos(i)
            return exp_pos

        if not (0 <= axis < self.num_ext_dim):
            raise ValueError(f"axis must be in [0, {self.num_ext_dim - 1}]")

        xp = self._engine.xp
        phys_axis = active_axes[axis]
        axes_to_trace = [0, 1, 2]
        axes_to_trace.pop(phys_axis)
        psi_sq_amp = xp.abs(self._amp) ** 2
        traced_out_psi = xp.sum(psi_sq_amp, axis=tuple(axes_to_trace))
        exp_pos = xp.sum(traced_out_psi * self._r[phys_axis][:, xp.newaxis])
        exp_pos *= self._volume_element()
        return exp_pos

    def var_pos(self, axis: int | None = None) -> NDArray[np.floating[Any]]:
        """Calculate position variance for one axis or all external axes."""
        active_axes = self._active_axis_indices()
        if axis is None:
            var_pos = np.empty((self.num_ext_dim,))
            for i in range(self.num_ext_dim):
                var_pos[i] = self.var_pos(i)
            return var_pos

        if not (0 <= axis < self.num_ext_dim):
            raise ValueError(f"axis must be in [0, {self.num_ext_dim - 1}]")

        xp = self._engine.xp
        phys_axis = active_axes[axis]
        axes_to_trace = [0, 1, 2]
        axes_to_trace.pop(phys_axis)
        psi_sq_amp = xp.abs(self._amp) ** 2
        traced_out_psi = xp.sum(psi_sq_amp, axis=tuple(axes_to_trace))
        variance_axis = (self._r[phys_axis] - self.exp_pos(axis)) ** 2
        var_pos = xp.sum(traced_out_psi * variance_axis[:, xp.newaxis])
        var_pos *= self._volume_element()
        return var_pos

    def normalize_to(self, n_const: float) -> None:
        """Normalize total wavefunction occupation to ``n_const``."""
        xp = self._engine.xp
        s = self._engine.inner_product(
            self._amp,
            self._amp,
            volume_element=self._volume_element(),
        )
        self._amp *= xp.sqrt(n_const / s)

    def state_occupation(
        self, nth_state: int | None = None
    ) -> NDArray[np.floating[Any]]:
        """Return occupation number of one internal state or all states."""
        if nth_state is None:
            occ = np.empty((self.num_int_dim,))
            for i in range(self.num_int_dim):
                occ[i] = self.state_occupation(i)
            return occ

        if not (0 <= nth_state < self.num_int_dim):
            raise ValueError(f"nth_state must be in [0, {self.num_int_dim - 1}]")

        return self._engine.state_occupation(
            self._amp,
            state_index=nth_state,
            volume_element=self._volume_element(),
        )

    def freely_propagate(
        self,
        steps: int,
        dt: float,
        options: PropagationOptions | None = None,
    ) -> None:
        """Propagate the wavefunction in time with ``V = 0``."""
        pytalises.propagator.freely_propagate(
            self,
            steps,
            dt,
            options=options,
        )

    def propagate(
        self,
        potential: BasePotential,
        steps: int,
        dt: float,
        *,
        variables: dict[str, Any] | None = None,
        options: PropagationOptions | None = None,
    ) -> None:
        """Propagate the wavefunction in time with a structured potential."""
        pytalises.propagator.propagate(
            self,
            potential,
            steps,
            dt,
            variables=variables,
            options=options,
        )
