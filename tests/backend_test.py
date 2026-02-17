import pytalises as pt


def test_default_backend_is_numpy():
    backend = pt.get_backend()
    assert backend.name == "numpy"


def test_wavefunction_uses_backend_instance():
    backend = pt.get_backend("numpy", num_threads=2)
    grid = pt.Grid(shape=(32,), extent=((-1, 1),))
    psi = pt.Wavefunction("exp(-x**2)", grid=grid, backend=backend)
    assert psi._backend is backend
