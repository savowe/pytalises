"""Structured potential definitions for pyTALISES v2."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


@dataclass(frozen=True)
class PotentialSpec:
    """Internal normalized potential representation."""

    expressions: tuple[str, ...]
    diag: bool


class BasePotential:
    """Base class for public potential specifications."""

    def to_spec(self, num_states: int) -> PotentialSpec:  # pragma: no cover - interface
        raise NotImplementedError


@dataclass(frozen=True)
class DiagonalPotential(BasePotential):
    """Diagonal potential with one expression per internal state."""

    expressions: tuple[str, ...]

    def __init__(self, expressions: str | Sequence[str]):
        if isinstance(expressions, str):
            expressions = (expressions,)
        object.__setattr__(self, "expressions", tuple(expressions))

    def to_spec(self, num_states: int) -> PotentialSpec:
        if len(self.expressions) != num_states:
            raise ValueError(
                "DiagonalPotential must contain exactly one expression per internal state. "
                f"Expected {num_states}, got {len(self.expressions)}."
            )
        return PotentialSpec(expressions=self.expressions, diag=True)


@dataclass(frozen=True)
class HermitianPotential(BasePotential):
    """Hermitian potential specified through lower-triangular elements."""

    lower_triangular: tuple[str, ...]

    def __init__(self, lower_triangular: Sequence[str]):
        object.__setattr__(self, "lower_triangular", tuple(lower_triangular))

    @property
    def num_states(self) -> int:
        num_v = len(self.lower_triangular)
        n = 0.5 * (math.sqrt(8 * num_v + 1) - 1)
        if not float(n).is_integer():
            raise ValueError(
                "Number of lower-triangular elements is invalid for a Hermitian matrix."
            )
        return int(n)

    @classmethod
    def from_lower_triangular(
        cls, lower_triangular: Sequence[str], num_states: int | None = None
    ) -> "HermitianPotential":
        pot = cls(lower_triangular)
        if num_states is not None and pot.num_states != num_states:
            raise ValueError(
                f"Lower triangular data defines {pot.num_states} states, expected {num_states}."
            )
        return pot

    @classmethod
    def from_matrix(cls, matrix: Sequence[Sequence[str]]) -> "HermitianPotential":
        n = len(matrix)
        if n == 0:
            raise ValueError("Potential matrix must not be empty.")
        if any(len(row) != n for row in matrix):
            raise ValueError("Potential matrix must be square.")

        lower: list[str] = []
        for i in range(n):
            for j in range(i + 1):
                if i == j:
                    lower.append(matrix[i][j])
                    continue
                lij = matrix[i][j]
                uji = matrix[j][i]
                if lij is None and uji is None:
                    raise ValueError(
                        f"Off-diagonal element ({i},{j})/( {j},{i}) must be provided."
                    )
                if lij is None:
                    lower.append(uji)
                elif uji is None:
                    lower.append(lij)
                elif lij == uji:
                    lower.append(lij)
                else:
                    raise ValueError(
                        "from_matrix currently expects symmetric symbolic entries for off-diagonal terms. "
                        f"Got ({i},{j})='{lij}' and ({j},{i})='{uji}'."
                    )
        return cls(tuple(lower))

    def to_spec(self, num_states: int) -> PotentialSpec:
        if self.num_states != num_states:
            raise ValueError(
                f"HermitianPotential defines {self.num_states} states but wavefunction has {num_states}."
            )
        return PotentialSpec(expressions=self.lower_triangular, diag=False)


def zero_potential(num_states: int) -> DiagonalPotential:
    """Return a diagonal zero-potential for all internal states."""
    return DiagonalPotential(["0"] * num_states)
