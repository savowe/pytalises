"""Grid definition for pyTALISES v2."""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class Grid:
    """Explicit spatial grid definition.

    Parameters
    ----------
    shape:
        Number of points per external dimension (1D-3D).
    extent:
        Spatial bounds per external dimension as ``(min, max)`` tuples.
    """

    shape: tuple[int, ...]
    extent: tuple[tuple[float, float], ...]

    def __post_init__(self) -> None:
        if not (1 <= len(self.shape) <= 3):
            raise ValueError("Grid supports 1 to 3 external dimensions.")
        if len(self.extent) != len(self.shape):
            raise ValueError("Grid.extent must have same length as Grid.shape.")
        for n in self.shape:
            if not isinstance(n, int) or n <= 0:
                raise ValueError("Grid.shape must contain positive integers.")
        for lo, hi in self.extent:
            if hi < lo:
                raise ValueError("Each extent tuple must satisfy max >= min.")

    @property
    def ndim(self) -> int:
        """Number of external dimensions."""
        return len(self.shape)

    @property
    def padded_shape(self) -> tuple[int, int, int]:
        """Shape padded to internal 3D representation."""
        return self.shape + (1,) * (3 - self.ndim)

    @property
    def padded_extent(self) -> tuple[tuple[float, float], tuple[float, float], tuple[float, float]]:
        """Extent padded to internal 3D representation."""
        return self.extent + ((0.0, 0.0),) * (3 - self.ndim)
