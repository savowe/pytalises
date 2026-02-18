"""Factory for private propagation engines."""

from __future__ import annotations

from pytalises.engine.cpu import CpuEngine


def create_engine(backend):
    """Create engine implementation for a backend instance."""
    if backend.name == "numpy":
        return CpuEngine(backend)

    if backend.name == "cupy":
        from pytalises.engine.gpu import GpuEngine

        return GpuEngine(backend)

    raise ValueError(f"No engine registered for backend '{backend.name}'.")
