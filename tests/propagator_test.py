import numpy as np
import pytalises as pt


def g1(n=64, a=-2.0, b=2.0):
    return pt.Grid(shape=(n,), extent=((a, b),))


def test_potential_static_and_linear_detection():
    static_potential = pt.DiagonalPotential(["x", "y"])
    p = pt.propagator.Propagator.Potential.from_potential(static_potential, {}, 2)
    assert p.static is True
    assert p.linear is True

    time_dependent = pt.DiagonalPotential(["t"])
    p = pt.propagator.Propagator.Potential.from_potential(time_dependent, {}, 1)
    assert p.static is False

    nonlinear = pt.DiagonalPotential(["psi0"])
    p = pt.propagator.Propagator.Potential.from_potential(nonlinear, {}, 1)
    assert p.static is False
    assert p.linear is False


def test_potential_internal_dimension_resolution():
    two_state_diag = pt.DiagonalPotential(["x", "2*x"])
    p_diag = pt.propagator.Propagator.Potential.from_potential(two_state_diag, {}, 2)
    assert p_diag.num_int_dim == 2

    two_state_nondiag = pt.HermitianPotential.from_lower_triangular(["x", "0", "y"])
    p_non = pt.propagator.Propagator.Potential.from_potential(two_state_nondiag, {}, 2)
    assert p_non.num_int_dim == 2


def test_wavefunction_normalization_and_moments():
    psi = pt.Wavefunction("exp(-x**2)", g1(64, -2, 2), normalize_const=1.0)

    r = psi.r
    k = psi.k
    assert isinstance(r, tuple)
    assert isinstance(k, tuple)
    assert r[0].shape == (64,)
    assert k[0].shape == (64,)

    exp_all = psi.exp_pos()
    exp_axis = psi.exp_pos(axis=0)
    np.testing.assert_allclose(exp_all, np.array([exp_axis]))

    var_all = psi.var_pos()
    var_axis = psi.var_pos(axis=0)
    np.testing.assert_allclose(var_all, np.array([var_axis]))

    psi.normalize_to(3.0)
    np.testing.assert_allclose(np.sum(psi.state_occupation()), 3.0)


def test_time_dependent_diagonal_potential():
    psi = pt.Wavefunction(
        ["exp(-x**2)", "0"],
        g1(64, -3, 3),
        normalize_const=1.0,
    )
    psi.propagate(
        potential=pt.DiagonalPotential(["t*x**2", "t*x**2"]),
        steps=5,
        dt=0.01,
    )
    np.testing.assert_almost_equal(np.sum(psi.state_occupation()), 1.0, decimal=5)


def test_time_dependent_nondiagonal_potential():
    psi = pt.Wavefunction(
        ["exp(-x**2)", "0"],
        g1(32, -2, 2),
        normalize_const=1.0,
    )
    psi.propagate(
        potential=pt.HermitianPotential.from_lower_triangular(["0", "sin(t)", "0"]),
        steps=3,
        dt=0.01,
    )
    np.testing.assert_almost_equal(np.sum(psi.state_occupation()), 1.0, decimal=5)


def test_multidimensional_wavefunction_properties():
    psi_2d = pt.Wavefunction(
        "exp(-x**2 - y**2)",
        pt.Grid(shape=(32, 32), extent=((-2, 2), (-2, 2))),
        normalize_const=1.0,
    )
    r = psi_2d.r
    k = psi_2d.k
    assert isinstance(r, tuple)
    assert isinstance(k, tuple)
    assert len(r) == 2
    assert len(k) == 2
    assert r[0].shape == (32,)
    assert r[1].shape == (32,)
    assert k[0].shape == (32,)
    assert k[1].shape == (32,)

    assert psi_2d.amp.shape == (32, 32, 1)
