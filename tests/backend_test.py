import pytest

import pytalises as pt


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


@pytest.mark.skipif(not pt.has_cupy(), reason="CuPy backend not available")
def test_wavefunction_can_run_on_cupy_backend():
    cupy = pytest.importorskip("cupy")
    grid = pt.Grid(shape=(16,), extent=((-1, 1),))
    psi = pt.Wavefunction("exp(-x**2)", grid=grid, backend="cupy")
    assert isinstance(psi.amp, cupy.ndarray)
    psi.freely_propagate(steps=2, dt=0.01)
    occ = psi.state_occupation()
    assert occ.shape == (1,)
