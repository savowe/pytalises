"""Module containing functionality for wavefunction propagation."""

from __future__ import annotations

from collections import defaultdict
from typing import TYPE_CHECKING, Any
import time

import numexpr as ne

from pytalises.engine import ExpressionEvaluator
from pytalises.options import PropagationOptions
from pytalises.potentials import BasePotential, zero_potential

if TYPE_CHECKING:
    from pytalises.wavefunction import Wavefunction

import pytalises.wavefunction


def propagate(
    psi: Wavefunction,
    potential: BasePotential,
    steps: int,
    dt: float,
    *,
    variables: dict[str, Any] | None = None,
    options: PropagationOptions | None = None,
) -> None:
    """Propagate a wavefunction in time using Strang splitting."""
    U = Propagator(
        psi,
        potential,
        variables=variables,
        options=options,
    )
    U.kinetic_prop(dt / 2)
    U.potential_prop(dt)
    for _ in range(steps - 1):
        U.kinetic_prop(dt)
        U.potential_prop(dt)
    U.kinetic_prop(dt / 2)


def freely_propagate(
    psi: Wavefunction,
    steps: int,
    dt: float,
    *,
    options: PropagationOptions | None = None,
) -> None:
    """Propagate a wavefunction in time with ``V = 0``."""
    U = Propagator(
        psi,
        potential=zero_potential(psi.num_int_dim),
        variables={},
        options=options,
    )
    for _ in range(steps):
        U.kinetic_prop(dt)


class Propagator:
    """Class for propagating instances of :class:`pytalises.Wavefunction`."""

    def __init__(
        self,
        psi: Wavefunction,
        potential: BasePotential,
        *,
        variables: dict[str, Any] | None = None,
        options: PropagationOptions | None = None,
    ):
        """Initialize propagator."""
        self.psi = psi
        self.options = options or PropagationOptions()

        if self.options.backend not in ("auto", self.psi._backend.name):
            raise ValueError(
                "PropagationOptions.backend must match Wavefunction backend. "
                f"Got options.backend='{self.options.backend}', "
                f"wavefunction backend='{self.psi._backend.name}'."
            )

        self._backend = self.psi._backend
        self._engine = self.psi._engine
        self._engine.set_num_threads(self.options.threads)
        self._evaluator = ExpressionEvaluator(backend_name=self._backend.name)

        self._profile_stages = bool(self.options.profile_stages)
        self._stage_seconds: dict[str, float] = defaultdict(float)
        self._stage_calls: dict[str, int] = defaultdict(int)

        self.v = self.Potential.from_potential(
            potential=potential,
            variables=variables or {},
            num_int_dim=psi.num_int_dim,
        )
        self._use_analytic_2x2 = (
            self.options.coupled_2x2_mode == "auto"
            and self.v.diag is False
            and self.psi.num_int_dim == 2
        )

        self._use_affine_precompute_diag = False
        self._use_affine_precompute_full = False
        self._affine_diag_a: Any | None = None
        self._affine_diag_b: Any | None = None
        self._affine_matrix_a: Any | None = None
        self._affine_matrix_b: Any | None = None

        assert isinstance(psi, pytalises.wavefunction.Wavefunction)
        assert self.v.num_int_dim == self.psi.num_int_dim
        assert self.psi._amp.shape[-1] == self.psi.num_int_dim

        self.V_eval_array = self._engine.zeros(
            psi.number_of_grid_points + (psi.num_int_dim, psi.num_int_dim),
            dtype="complex128",
        )
        self.V_eval_diag_array = self._engine.zeros(
            psi.number_of_grid_points + (psi.num_int_dim,),
            dtype="complex128",
        )
        self.V_eval_eigval_array = self._engine.zeros(
            psi.number_of_grid_points + (psi.num_int_dim,),
            dtype="complex128",
        )
        self.V_eval_eigvec_array = self._engine.zeros(
            psi.number_of_grid_points + (psi.num_int_dim, psi.num_int_dim),
            dtype="complex128",
        )

        self.psi.construct_FFT(self.options.threads, self.options.fftw_flags)

        if self.v.diag is True:
            self.prop_method = self.diag_potential_prop
        else:
            self.prop_method = self.nondiag_potential_prop

        self._configure_affine_time_precompute()

        if self.v.static is True:
            if self.v.diag is True:
                self.eval_diag_V()
            else:
                self.eval_V()
                if not self._use_analytic_2x2:
                    self._refresh_eigendecomposition()

    def _to_scalar_float(self, value: Any) -> float:
        if hasattr(value, "item"):
            return float(value.item())
        return float(value)

    def _max_abs(self, value: Any) -> float:
        xp = self._engine.xp
        return self._to_scalar_float(xp.max(xp.abs(value)))

    def _try_affine_decompose_time(
        self,
        expression: str,
        local_scope: dict[str, Any],
    ) -> tuple[Any, Any] | None:
        eval_t0 = self._evaluator.eval(
            expression,
            local_dict=local_scope,
            global_dict={"t": 0.0},
        )
        eval_t1 = self._evaluator.eval(
            expression,
            local_dict=local_scope,
            global_dict={"t": 1.0},
        )
        eval_t2 = self._evaluator.eval(
            expression,
            local_dict=local_scope,
            global_dict={"t": 2.0},
        )

        a = eval_t0
        b = eval_t1 - eval_t0
        residual = eval_t2 - (a + 2.0 * b)

        scale = max(
            1.0,
            self._max_abs(eval_t0),
            self._max_abs(eval_t1),
            self._max_abs(eval_t2),
        )
        error = self._max_abs(residual)

        if error <= 1e-12 + 1e-10 * scale:
            return a, b
        return None

    def _build_affine_diag_precompute(self) -> bool:
        local_scope, _ = self._evaluation_scope()
        a_diag = self._engine.zeros(self.V_eval_diag_array.shape, dtype="complex128")
        b_diag = self._engine.zeros(self.V_eval_diag_array.shape, dtype="complex128")

        for i, expression in enumerate(self.v.potential_strings):
            coeffs = self._try_affine_decompose_time(expression, local_scope)
            if coeffs is None:
                return False
            a_diag[:, :, :, i] = coeffs[0]
            b_diag[:, :, :, i] = coeffs[1]

        self._affine_diag_a = a_diag
        self._affine_diag_b = b_diag
        self._use_affine_precompute_diag = True
        return True

    def _build_affine_full_precompute(self) -> bool:
        local_scope, _ = self._evaluation_scope()
        a_full = self._engine.zeros(self.V_eval_array.shape, dtype="complex128")
        b_full = self._engine.zeros(self.V_eval_array.shape, dtype="complex128")

        k = 0
        for i in range(self.psi.num_int_dim):
            for j in range(i, self.psi.num_int_dim):
                coeffs = self._try_affine_decompose_time(
                    self.v.potential_strings[k],
                    local_scope,
                )
                if coeffs is None:
                    return False

                a_ji, b_ji = coeffs
                a_full[:, :, :, j, i] = a_ji
                b_full[:, :, :, j, i] = b_ji
                if i != j:
                    a_full[:, :, :, i, j] = self._engine.xp.conjugate(a_ji)
                    b_full[:, :, :, i, j] = self._engine.xp.conjugate(b_ji)
                k += 1

        self._affine_matrix_a = a_full
        self._affine_matrix_b = b_full
        self._use_affine_precompute_full = True
        return True

    def _configure_affine_time_precompute(self) -> None:
        if self.options.potential_precompute_mode != "auto":
            return
        if self.v.static is True:
            return
        if self.v.linear is False:
            return

        if self.v.diag is True:
            self._build_affine_diag_precompute()
        else:
            self._build_affine_full_precompute()

    def _record_stage(self, name: str, elapsed: float) -> None:
        if not self._profile_stages:
            return
        self._stage_seconds[name] += elapsed
        self._stage_calls[name] += 1

    def _start_stage_timer(self) -> float:
        if not self._profile_stages:
            return 0.0
        self._engine.synchronize()
        return time.perf_counter()

    def _stop_stage_timer(self, name: str, start: float) -> None:
        if not self._profile_stages:
            return
        self._engine.synchronize()
        self._record_stage(name, time.perf_counter() - start)

    def stage_timings(self) -> dict[str, dict[str, float | int]]:
        """Return accumulated stage timings for this propagator instance."""
        if not self._profile_stages:
            return {}
        return {
            name: {
                "seconds": float(self._stage_seconds[name]),
                "calls": int(self._stage_calls[name]),
            }
            for name in sorted(self._stage_seconds)
        }

    def potential_prop(self, dt: float) -> None:
        """Apply potential propagator step."""
        self.prop_method(dt)

    def _refresh_eigendecomposition(self) -> None:
        t0 = self._start_stage_timer()
        eigvals, eigvecs = self._engine.eigendecompose_hermitian(self.V_eval_array)
        self.V_eval_eigval_array[...] = eigvals
        self.V_eval_eigvec_array[...] = eigvecs
        self._stop_stage_timer("_refresh_eigendecomposition", t0)

    def nondiag_potential_prop(self, dt: float) -> None:
        """Apply potential step for non-diagonal potentials."""
        if self.v.static is False:
            self.eval_V()
            if not self._use_analytic_2x2:
                self._refresh_eigendecomposition()

        if self._use_analytic_2x2:
            t0 = self._start_stage_timer()
            self._engine.apply_coupled_phase_2x2(
                self.psi._amp,
                matrix=self.V_eval_array,
                dt=dt,
            )
            self._stop_stage_timer("apply_coupled_phase_analytic_2x2", t0)
            return

        t0 = self._start_stage_timer()
        self._engine.apply_coupled_phase(
            self.psi._amp,
            eigvals=self.V_eval_eigval_array,
            eigvecs=self.V_eval_eigvec_array,
            dt=dt,
        )
        self._stop_stage_timer("apply_coupled_phase", t0)

    def diag_potential_prop(self, dt: float) -> None:
        """Apply potential step for diagonal potentials."""
        if self.v.static is False:
            self.eval_diag_V()
        t0 = self._start_stage_timer()
        self._engine.apply_diagonal_phase(
            self.psi._amp,
            diagonal=self.V_eval_diag_array,
            dt=dt,
        )
        self._stop_stage_timer("apply_diagonal_phase", t0)

    def kinetic_prop(self, dt: float) -> None:
        """Perform kinetic propagation step in reciprocal space."""
        t0 = self._start_stage_timer()
        self.psi.fft()
        self._engine.apply_kinetic_phase(
            self.psi._amp,
            kmesh=(self.psi.kmesh[0], self.psi.kmesh[1], self.psi.kmesh[2]),
            alpha=self.psi.alpha,
            dt=dt,
        )
        self.psi.ifft()
        self.psi.t += dt
        self._stop_stage_timer("kinetic_prop", t0)

    def _evaluation_scope(self) -> tuple[dict[str, Any], dict[str, Any]]:
        local_scope = {**self.v.variables, **self.psi.default_var_dict}
        global_scope = {"t": self.psi.t}
        return local_scope, global_scope

    def eval_V(self) -> None:
        """Evaluate full potential matrix on the complete spatial grid."""
        t0 = self._start_stage_timer()

        if self._use_affine_precompute_full:
            assert self._affine_matrix_a is not None
            assert self._affine_matrix_b is not None
            self.V_eval_array[...] = self._affine_matrix_b
            self.V_eval_array *= self.psi.t
            self.V_eval_array += self._affine_matrix_a
            self._stop_stage_timer("eval_V", t0)
            return

        self.V_eval_array[...] = 0
        local_scope, global_scope = self._evaluation_scope()
        k = 0
        for i in range(self.psi.num_int_dim):
            for j in range(i, self.psi.num_int_dim):
                eval_ji = self._evaluator.eval(
                    self.v.potential_strings[k],
                    local_dict=local_scope,
                    global_dict=global_scope,
                )
                self.V_eval_array[:, :, :, j, i] = eval_ji
                if i != j:
                    self.V_eval_array[:, :, :, i, j] = self._engine.xp.conjugate(eval_ji)
                k += 1
        self._stop_stage_timer("eval_V", t0)

    def eval_diag_V(self) -> None:
        """Evaluate diagonal potential matrix elements."""
        t0 = self._start_stage_timer()

        if self._use_affine_precompute_diag:
            assert self._affine_diag_a is not None
            assert self._affine_diag_b is not None
            self.V_eval_diag_array[...] = self._affine_diag_b
            self.V_eval_diag_array *= self.psi.t
            self.V_eval_diag_array += self._affine_diag_a
            self._stop_stage_timer("eval_diag_V", t0)
            return

        local_scope, global_scope = self._evaluation_scope()
        for i in range(self.psi.num_int_dim):
            self.V_eval_diag_array[:, :, :, i] = self._evaluator.eval(
                self.v.potential_strings[i],
                local_dict=local_scope,
                global_dict=global_scope,
            )
        self._stop_stage_timer("eval_diag_V", t0)

    class Potential:
        """Simple container for potential metadata."""

        @classmethod
        def from_potential(
            cls,
            potential: BasePotential,
            variables: dict[str, Any],
            num_int_dim: int,
        ) -> "Propagator.Potential":
            if not isinstance(potential, BasePotential):
                raise TypeError(
                    "potential must be a structured potential object "
                    "(e.g. DiagonalPotential or HermitianPotential)."
                )
            spec = potential.to_spec(num_states=num_int_dim)
            return cls(
                potential_string=list(spec.expressions),
                variables=variables,
                diag=spec.diag,
            )

        def __init__(
            self,
            potential_string: list[str],
            variables: dict[str, Any] | None = None,
            diag: bool = False,
        ) -> None:
            if variables is None:
                variables = {}

            self.potential_strings = potential_string
            self.num_v = len(self.potential_strings)
            self.variables = variables

            for pot_string in self.potential_strings:
                potential_nex = ne.NumExpr(pot_string)
                try:
                    potential_nex.input_names.index("t")
                    self.static = False
                except ValueError:
                    self.static = True
                if self.static is False:
                    break

            for i in range(len(self.potential_strings)):
                for pot_string in self.potential_strings:
                    potential_nex = ne.NumExpr(pot_string)
                    try:
                        potential_nex.input_names.index("psi" + str(i))
                        self.linear = False
                    except ValueError:
                        self.linear = True
                    if self.linear is False:
                        self.static = False
                        break
                if self.linear is False:
                    self.static = False
                    break

            self.diag = diag
            if diag is False:
                self.num_int_dim = 1 / 2 * (self._sqrt(8 * self.num_v + 1) - 1)
                assert (
                    self.num_int_dim.is_integer()
                ), "Number of potential matrix elements incorrect"
                self.num_int_dim = int(self.num_int_dim)
            else:
                self.num_int_dim = len(self.potential_strings)

        @staticmethod
        def _sqrt(x: float) -> float:
            return float(x) ** 0.5
