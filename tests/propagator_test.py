import pytalises as pt
import numpy as np


def test_potential_static_and_linear_detection():
    static_potential = pt.propagator.Propagator.Potential(["x", "y"], diag=True)
    assert static_potential.static is True
    assert static_potential.linear is True

    time_dependent = pt.propagator.Propagator.Potential(["t"], diag=True)
    assert time_dependent.static is False

    nonlinear = pt.propagator.Propagator.Potential(["psi0"], diag=True)
    assert nonlinear.static is False
    assert nonlinear.linear is False


def test_potential_internal_dimension_resolution():
    two_state_diag = pt.propagator.Propagator.Potential(["x", "2*x"], diag=True)
    assert two_state_diag.num_int_dim == 2

    two_state_nondiag = pt.propagator.Propagator.Potential(["x", "0", "y"], diag=False)
    assert two_state_nondiag.num_int_dim == 2


def test_wavefunction_normalization_and_moments():
    psi = pt.Wavefunction("exp(-x**2)", (64,), (-2, 2), normalize_const=1.0)

    exp_all = psi.exp_pos()
    exp_axis = psi.exp_pos(axis=0)
    np.testing.assert_allclose(exp_all, np.array([exp_axis]))

    var_all = psi.var_pos()
    var_axis = psi.var_pos(axis=0)
    np.testing.assert_allclose(var_all, np.array([var_axis]))

    psi.normalize_to(3.0)
    np.testing.assert_allclose(np.sum(psi.state_occupation()), 3.0)
