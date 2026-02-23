"""Expression evaluation helpers for propagation internals.

This module intentionally keeps expression handling separate from backend kernels.
"""

from __future__ import annotations

from typing import Any

import numexpr as ne


class ExpressionEvaluator:
    """Evaluate symbolic expressions on CPU/GPU array backends."""

    def __init__(self, *, backend_name: str):
        self.backend_name = backend_name
        self._gpu_scope = None

        if backend_name == "cupy":
            import cupy as cp  # optional dependency

            self._gpu_scope = {
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

    def eval(
        self,
        expr: str,
        *,
        local_dict: dict[str, Any],
        global_dict: dict[str, Any] | None = None,
    ) -> Any:
        """Evaluate expression against backend-native arrays."""
        if self.backend_name == "numpy":
            return ne.evaluate(
                expr,
                local_dict=local_dict,
                global_dict=global_dict or {},
                order="C",
            )

        if self.backend_name == "cupy":
            scope = dict(self._gpu_scope or {})
            if global_dict:
                scope.update(global_dict)
            scope.update(local_dict)

            # CuPy ufunc dispatch may import internals dynamically.
            safe_builtins = {"__import__": __import__}
            return eval(expr, {"__builtins__": safe_builtins}, scope)

        raise ValueError(f"Unsupported backend for expression evaluation: {self.backend_name}")
