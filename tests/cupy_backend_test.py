import numpy as np
import pytest

import pytalises as pt


pytestmark = pytest.mark.skipif(not pt.has_cupy(), reason="CuPy backend not available")


def g1(n=64, a=-2.0, b=2.0):
    return pt.Grid(shape=(n,), extent=((a, b),))


def _to_numpy(x):
    import cupy as cp

    if isinstance(x, cp.ndarray):
        return cp.asnumpy(x)
    return np.asarray(x)


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


def test_cupy_matches_numpy_reference_dynamics():
    grid = g1(96, -4, 4)
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*sin(t)",
            "0.1*x**2 - 0.05*t",
        ]
    )

    psi_np = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="numpy")
    psi_cp = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="cupy")

    steps = 10
    dt = 0.005
    psi_np.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(backend="numpy"),
    )
    psi_cp.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(backend="cupy"),
    )

    amp_np = _to_numpy(psi_np.amp)
    amp_cp = _to_numpy(psi_cp.amp)

    np.testing.assert_allclose(
        np.abs(amp_cp) ** 2,
        np.abs(amp_np) ** 2,
        rtol=5e-5,
        atol=5e-7,
    )
    np.testing.assert_allclose(
        _to_numpy(psi_cp.state_occupation()),
        _to_numpy(psi_np.state_occupation()),
        rtol=5e-5,
        atol=5e-7,
    )


def test_cupy_backend_einsum_out_semantics():
    cp = pytest.importorskip("cupy")
    backend = pt.get_backend("cupy")

    a = cp.arange(6, dtype=cp.float64).reshape(2, 3)
    b = cp.arange(3, dtype=cp.float64)
    out = cp.empty((2,), dtype=cp.float64)

    result = backend.einsum("ij,j->i", a, b, out=out)

    np.testing.assert_allclose(cp.asnumpy(out), cp.asnumpy(cp.einsum("ij,j->i", a, b)))
    assert result is out


def test_cupy_analytic_2x2_matches_eigh_reference():
    grid = g1(96, -4, 4)
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*sin(t)",
            "0.1*x**2 - 0.05*t",
        ]
    )

    psi_auto = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="cupy")
    psi_eigh = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="cupy")

    steps = 10
    dt = 0.005
    psi_auto.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(backend="cupy", coupled_2x2_mode="auto"),
    )
    psi_eigh.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(backend="cupy", coupled_2x2_mode="eigh"),
    )

    amp_auto = _to_numpy(psi_auto.amp)
    amp_eigh = _to_numpy(psi_eigh.amp)

    np.testing.assert_allclose(
        np.abs(amp_auto) ** 2,
        np.abs(amp_eigh) ** 2,
        rtol=5e-5,
        atol=5e-7,
    )
    np.testing.assert_allclose(
        _to_numpy(psi_auto.state_occupation()),
        _to_numpy(psi_eigh.state_occupation()),
        rtol=5e-5,
        atol=5e-7,
    )


def test_cupy_complex64_mode_matches_complex128_reference():
    grid = g1(96, -4, 4)
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*sin(t)",
            "0.1*x**2 - 0.05*t",
        ]
    )

    psi64 = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="cupy")
    psi128 = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="cupy")

    steps = 10
    dt = 0.005

    psi64.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(
            backend="cupy",
            dtype="complex64",
            coupled_2x2_mode="auto",
        ),
    )
    psi128.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(
            backend="cupy",
            dtype="complex128",
            coupled_2x2_mode="auto",
        ),
    )

    amp64 = _to_numpy(psi64.amp)
    amp128 = _to_numpy(psi128.amp)

    assert amp64.dtype == np.complex64

    np.testing.assert_allclose(
        np.abs(amp64) ** 2,
        np.abs(amp128) ** 2,
        rtol=2e-3,
        atol=2e-5,
    )
    np.testing.assert_allclose(
        _to_numpy(psi64.state_occupation()),
        _to_numpy(psi128.state_occupation()),
        rtol=2e-3,
        atol=2e-5,
    )
