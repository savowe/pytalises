import numpy as np
import pytest

import pytalises as pt


def _to_numpy(x):
    if hasattr(x, "get"):
        return np.asarray(x.get())
    return np.asarray(x)


def _backend_runtime_available(name: str) -> bool:
    if name == "cupy":
        return pt.has_cupy()
    return True


def _run_reference_case(backend: str):
    grid = pt.Grid(shape=(96,), extent=((-4.0, 4.0),))
    initial = ["exp(-x**2)", "0.2*exp(-(x-0.5)**2)"]
    potential = pt.HermitianPotential.from_lower_triangular(
        [
            "0.2*x**2 + 0.1*t",
            "0.03*cos(x) + 0.02*sin(t)",
            "0.1*x**2 - 0.05*t",
        ]
    )

    psi = pt.Wavefunction(initial, grid, normalize_const=1.0, backend=backend)
    psi.propagate(
        potential=potential,
        steps=10,
        dt=0.005,
        options=pt.PropagationOptions(backend=backend),
    )
    return psi


CANDIDATE_BACKENDS = [name for name in pt.available_backends() if name != "numpy"]
if not CANDIDATE_BACKENDS:
    CANDIDATE_BACKENDS = [
        pytest.param(
            "_none_",
            marks=pytest.mark.skip(reason="No non-NumPy backends are importable"),
        )
    ]


@pytest.mark.parametrize("backend_name", CANDIDATE_BACKENDS)
def test_backend_matches_numpy_reference(backend_name: str):
    if not _backend_runtime_available(backend_name):
        pytest.skip(f"{backend_name} runtime is not available")

    psi_np = _run_reference_case("numpy")
    psi_other = _run_reference_case(backend_name)

    amp_np = _to_numpy(psi_np.amp)
    amp_other = _to_numpy(psi_other.amp)

    np.testing.assert_allclose(
        np.abs(amp_other) ** 2,
        np.abs(amp_np) ** 2,
        rtol=5e-5,
        atol=5e-7,
    )
    np.testing.assert_allclose(
        _to_numpy(psi_other.state_occupation()),
        _to_numpy(psi_np.state_occupation()),
        rtol=5e-5,
        atol=5e-7,
    )


def test_numpy_backend_is_available():
    assert "numpy" in pt.available_backends()
