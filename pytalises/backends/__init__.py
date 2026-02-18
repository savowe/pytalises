"""Backend registry and auto-selection."""

from __future__ import annotations

from typing import Optional, Union

from .base import Backend
from .numpy_backend import NumpyBackend

_BACKENDS = {
    "numpy": NumpyBackend,
}

try:  # optional dependency
    from .cupy_backend import CupyBackend
except ImportError:  # pragma: no cover - depends on optional runtime
    CupyBackend = None
else:
    _BACKENDS["cupy"] = CupyBackend

_default_backend: Optional[Backend] = None


def register_backend(name: str, backend_class):
    """Register a backend class under ``name``."""
    _BACKENDS[name] = backend_class


def available_backends() -> tuple[str, ...]:
    """Return names of currently importable backends."""
    return tuple(sorted(_BACKENDS.keys()))


def has_cupy() -> bool:
    """Return ``True`` if a usable CuPy backend is available."""
    if "cupy" not in _BACKENDS:
        return False
    try:
        import cupy as cp

        return cp.cuda.runtime.getDeviceCount() > 0
    except Exception:
        return False


def _auto_detect_backend() -> str:
    if has_cupy():
        return "cupy"
    return "numpy"


def get_backend(name: Union[str, Backend, None] = None, **kwargs) -> Backend:
    """Return backend instance by name or pass-through existing instance."""
    global _default_backend

    if isinstance(name, Backend):
        return name

    if name is None or name == "auto":
        if _default_backend is not None and name is None:
            return _default_backend
        name = _auto_detect_backend()

    if name not in _BACKENDS:
        available = sorted(_BACKENDS.keys())
        raise ValueError(f"Unknown backend '{name}'. Available: {available}")

    if name == "cupy" and not has_cupy():
        raise ValueError(
            "CuPy backend requested but no usable CUDA device is available. "
            "Install CuPy with CUDA support and ensure at least one visible GPU."
        )

    return _BACKENDS[name](**kwargs)


def set_default_backend(backend: Union[str, Backend], **kwargs) -> None:
    """Set global default backend instance."""
    global _default_backend
    _default_backend = get_backend(backend, **kwargs)
