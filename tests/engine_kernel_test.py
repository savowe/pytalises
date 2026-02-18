import numpy as np

import pytalises as pt
from pytalises.engine import create_engine


def test_coupled_kernel_matches_manual_matrix_update():
    backend = pt.get_backend("numpy")
    engine = create_engine(backend)

    amp = backend.asarray([[[[1.0 + 0.0j, 0.25 - 0.1j]]]], dtype="complex128")
    eigvals = backend.asarray([[[[0.3, -0.2]]]], dtype="complex128")

    s = 1 / np.sqrt(2)
    eigvec = np.array([[s, s], [s, -s]], dtype=np.complex128)
    eigvecs = backend.asarray([[[eigvec]]], dtype="complex128")

    dt = 0.7

    expected = np.asarray(amp).copy()
    propagator = eigvec @ np.diag(np.exp(-1j * np.array([0.3, -0.2]) * dt)) @ np.conjugate(
        eigvec.T
    )
    expected[0, 0, 0, :] = propagator @ expected[0, 0, 0, :]

    engine.apply_coupled_phase(amp, eigvals=eigvals, eigvecs=eigvecs, dt=dt)

    np.testing.assert_allclose(np.asarray(amp), expected, atol=1e-12, rtol=1e-12)
