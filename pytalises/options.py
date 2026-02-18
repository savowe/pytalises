"""Runtime propagation options for pyTALISES v2."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PropagationOptions:
    """Runtime configuration for propagation.

    Notes
    -----
    ``backend`` and ``dtype`` are part of the v2 public surface and may be
    ignored by the current NumPy backend implementation until additional
    backends are integrated.
    """

    backend: str = "auto"
    dtype: str = "complex128"
    threads: int = 1
    fftw_flags: tuple[str, ...] = (
        "FFTW_ESTIMATE",
        "FFTW_DESTROY_INPUT",
    )
    profile_stages: bool = False
    coupled_2x2_mode: str = "auto"
    potential_precompute_mode: str = "off"

    def __post_init__(self) -> None:
        if self.threads < 1:
            raise ValueError("PropagationOptions.threads must be >= 1")
        if self.coupled_2x2_mode not in {"auto", "eigh"}:
            raise ValueError(
                "PropagationOptions.coupled_2x2_mode must be one of: 'auto', 'eigh'"
            )
        if self.potential_precompute_mode not in {"off", "auto"}:
            raise ValueError(
                "PropagationOptions.potential_precompute_mode must be one of: 'off', 'auto'"
            )
