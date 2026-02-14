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

    # Test r and k properties for 1D wavefunction
    r = psi.r
    k = psi.k
    assert isinstance(r, np.ndarray)
    assert isinstance(k, np.ndarray)
    assert r.shape == (64,)
    assert k.shape == (64,)

    exp_all = psi.exp_pos()
    exp_axis = psi.exp_pos(axis=0)
    np.testing.assert_allclose(exp_all, np.array([exp_axis]))

    var_all = psi.var_pos()
    var_axis = psi.var_pos(axis=0)
    np.testing.assert_allclose(var_all, np.array([var_axis]))

    psi.normalize_to(3.0)
    np.testing.assert_allclose(np.sum(psi.state_occupation()), 3.0)


def test_time_dependent_diagonal_potential():
    """Test propagation with a time-dependent diagonal potential."""
    psi = pt.Wavefunction(
        ["exp(-x**2)", "0"],
        (64,),
        (-3, 3),
        normalize_const=1.0,
    )
    # Time-dependent potential V = t * x^2 (harmonic with time-varying strength)
    psi.propagate(
        potential=["t*x**2", "t*x**2"],
        num_time_steps=5,
        delta_t=0.01,
        diag=True,
    )
    # Check normalization is preserved
    np.testing.assert_almost_equal(np.sum(psi.state_occupation()), 1.0, decimal=5)


def test_time_dependent_nondiagonal_potential():
    """Test propagation with a time-dependent nondiagonal potential."""
    psi = pt.Wavefunction(
        ["exp(-x**2)", "0"],
        (32,),
        (-2, 2),
        normalize_const=1.0,
    )
    # Time-dependent coupling: V01 = sin(t)
    # Lower triangular: [V00, V10, V11]
    psi.propagate(
        potential=["0", "sin(t)", "0"],
        num_time_steps=3,
        delta_t=0.01,
        diag=False,
    )
    # Normalization should be preserved
    np.testing.assert_almost_equal(np.sum(psi.state_occupation()), 1.0, decimal=5)


def test_multidimensional_wavefunction_properties():
    """Test r, k properties for multi-dimensional wavefunctions."""
    # 2D wavefunction
    psi_2d = pt.Wavefunction(
        "exp(-x**2 - y**2)",
        (32, 32),
        [(-2, 2), (-2, 2)],
        normalize_const=1.0,
    )
    # r and k should return lists for multi-dim
    r = psi_2d.r
    k = psi_2d.k
    assert isinstance(r, list)
    assert isinstance(k, list)
    assert len(r) == 3  # x, y, z (z is trivial)
    assert len(k) == 3
    # Verify the arrays have the right shapes
    assert r[0].shape == (32,)
    assert r[1].shape == (32,)
    assert k[0].shape == (32,)
    assert k[1].shape == (32,)

    # amp should be squeezed
    assert psi_2d.amp.shape == (32, 32)  # squeezed from (32, 32, 1, 1)
