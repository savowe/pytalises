import warnings

import numpy as np
import pytalises as pt


def test_legacy_wavefunction_and_propagation():
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        psi = pt.legacy.Wavefunction("exp(-x**2)", (64,), (-2, 2), normalize_const=1.0)
        psi.freely_propagate(num_time_steps=1, delta_t=0.01)
        psi.propagate("x**2", num_time_steps=1, delta_t=0.01, diag=True)

    assert any(issubclass(item.category, DeprecationWarning) for item in w)
    np.testing.assert_almost_equal(np.sum(psi.state_occupation()), 1.0, decimal=5)


def test_legacy_module_level_wrappers():
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", DeprecationWarning)
        psi = pt.legacy.Wavefunction(
            ["exp(-x**2)", "0"], (64,), (-2, 2), normalize_const=1.0
        )
        pt.legacy.propagate(
            psi,
            ["0", "sin(t)", "0"],
            num_time_steps=2,
            delta_t=0.01,
            diag=False,
        )
    np.testing.assert_almost_equal(np.sum(psi.state_occupation()), 1.0, decimal=5)
