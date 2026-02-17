import numpy as np
import pytest

import pytalises as pt


pytestmark = pytest.mark.skipif(not pt.has_cupy(), reason="CuPy backend not available")


def g1(n=64, a=-2.0, b=2.0):
    return pt.Grid(shape=(n,), extent=((a, b),))


def test_cupy_free_propagation_conserves_norm():
    psi = pt.Wavefunction(
        ["exp(-x**2)"],
        g1(64, -3, 3),
        normalize_const=1.0,
        backend="cupy",
    )
    psi.freely_propagate(
        steps=5,
        dt=0.01,
        options=pt.PropagationOptions(backend="cupy"),
    )
    np.testing.assert_allclose(np.sum(psi.state_occupation()), 1.0, atol=1e-5)


def test_cupy_diagonal_time_dependent_potential():
    psi = pt.Wavefunction(
        ["exp(-x**2)", "0"],
        g1(64, -3, 3),
        normalize_const=1.0,
        backend="cupy",
    )
    psi.propagate(
        potential=pt.DiagonalPotential(["t*x**2", "t*x**2"]),
        steps=4,
        dt=0.01,
        options=pt.PropagationOptions(backend="cupy"),
    )
    np.testing.assert_allclose(np.sum(psi.state_occupation()), 1.0, atol=1e-5)


def test_cupy_nondiagonal_potential():
    psi = pt.Wavefunction(
        ["exp(-x**2)", "0"],
        g1(32, -2, 2),
        normalize_const=1.0,
        backend="cupy",
    )
    psi.propagate(
        potential=pt.HermitianPotential.from_lower_triangular(["0", "sin(t)", "0"]),
        steps=3,
        dt=0.01,
        options=pt.PropagationOptions(backend="cupy"),
    )
    np.testing.assert_allclose(np.sum(psi.state_occupation()), 1.0, atol=1e-5)
