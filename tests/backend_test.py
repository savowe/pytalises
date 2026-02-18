import numpy as np
import pytest

import pytalises as pt
import pytalises.backends as backends


def test_default_backend_is_registered_backend():
    backend = pt.get_backend()
    assert backend.name in {"numpy", "cupy"}


def test_wavefunction_uses_backend_instance():
    backend = pt.get_backend("numpy", num_threads=2)
    grid = pt.Grid(shape=(32,), extent=((-1, 1),))
    psi = pt.Wavefunction("exp(-x**2)", grid=grid, backend=backend)
    assert psi._backend is backend


def test_requesting_cupy_backend_matches_availability():
    if pt.has_cupy():
        backend = pt.get_backend("cupy")
        assert backend.name == "cupy"
    else:
        with pytest.raises(ValueError):
            pt.get_backend("cupy")


@pytest.mark.skipif("cupy" not in pt.available_backends(), reason="CuPy backend module not importable")
def test_requesting_cupy_backend_fails_fast_without_visible_device(monkeypatch):
    monkeypatch.setattr(backends, "has_cupy", lambda: False)
    with pytest.raises(ValueError, match="no usable CUDA device"):
        backends.get_backend("cupy")


@pytest.mark.skipif(not pt.has_cupy(), reason="CuPy backend not available")
def test_wavefunction_can_run_on_cupy_backend():
    cupy = pytest.importorskip("cupy")
    grid = pt.Grid(shape=(16,), extent=((-1, 1),))
    psi = pt.Wavefunction("exp(-x**2)", grid=grid, backend="cupy")
    assert isinstance(psi.amp, cupy.ndarray)
    psi.freely_propagate(steps=2, dt=0.01)
    occ = psi.state_occupation()
    assert occ.shape == (1,)


def test_core_propagation_path_does_not_use_backend_einsum(monkeypatch):
    grid = pt.Grid(shape=(32,), extent=((-2, 2),))
    psi = pt.Wavefunction("exp(-x**2)", grid=grid, backend="numpy")

    def _unexpected_einsum(*args, **kwargs):
        raise AssertionError("Propagation should not call backend.einsum in v2 core path")

    monkeypatch.setattr(psi._backend, "einsum", _unexpected_einsum)

    psi.freely_propagate(steps=2, dt=0.01)
    psi.propagate(
        potential=pt.DiagonalPotential("0.1*x**2"),
        steps=2,
        dt=0.01,
        options=pt.PropagationOptions(backend="numpy"),
    )


def test_coupled_2x2_mode_rejects_invalid_value():
    with pytest.raises(ValueError, match="coupled_2x2_mode"):
        pt.PropagationOptions(coupled_2x2_mode="invalid")


def test_potential_precompute_mode_rejects_invalid_value():
    with pytest.raises(ValueError, match="potential_precompute_mode"):
        pt.PropagationOptions(potential_precompute_mode="invalid")


def test_analytic_2x2_path_matches_eigh_reference_numpy():
    grid = pt.Grid(shape=(64,), extent=((-3, 3),))
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*sin(t)",
            "0.1*x**2 - 0.05*t",
        ]
    )

    psi_auto = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="numpy")
    psi_eigh = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="numpy")

    steps = 12
    dt = 0.005

    psi_auto.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(backend="numpy", coupled_2x2_mode="auto"),
    )
    psi_eigh.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(backend="numpy", coupled_2x2_mode="eigh"),
    )

    np.testing.assert_allclose(
        np.abs(np.asarray(psi_auto.amp)) ** 2,
        np.abs(np.asarray(psi_eigh.amp)) ** 2,
        rtol=1e-10,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.asarray(psi_auto.state_occupation()),
        np.asarray(psi_eigh.state_occupation()),
        rtol=1e-10,
        atol=1e-12,
    )


def test_affine_potential_precompute_matches_off_reference_numpy():
    grid = pt.Grid(shape=(64,), extent=((-3, 3),))
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*t",
            "0.1*x**2 - 0.05*t",
        ]
    )

    psi_off = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="numpy")
    psi_auto = pt.Wavefunction(initial, grid, normalize_const=1.0, backend="numpy")

    steps = 12
    dt = 0.005

    psi_off.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(
            backend="numpy",
            coupled_2x2_mode="auto",
            potential_precompute_mode="off",
        ),
    )
    psi_auto.propagate(
        potential=potential,
        steps=steps,
        dt=dt,
        options=pt.PropagationOptions(
            backend="numpy",
            coupled_2x2_mode="auto",
            potential_precompute_mode="auto",
        ),
    )

    np.testing.assert_allclose(
        np.abs(np.asarray(psi_auto.amp)) ** 2,
        np.abs(np.asarray(psi_off.amp)) ** 2,
        rtol=1e-10,
        atol=1e-12,
    )
    np.testing.assert_allclose(
        np.asarray(psi_auto.state_occupation()),
        np.asarray(psi_off.state_occupation()),
        rtol=1e-10,
        atol=1e-12,
    )
