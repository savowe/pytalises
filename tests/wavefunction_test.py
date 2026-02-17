import itertools

import numpy as np
import pytalises as pt


def g1(n=64, a=-2.0, b=2.0):
    return pt.Grid(shape=(n,), extent=((a, b),))


def g2(nx=32, ny=32, a=-2.0, b=2.0):
    return pt.Grid(shape=(nx, ny), extent=((a, b), (a, b)))


def g3(n=16, a=-2.0, b=2.0):
    return pt.Grid(shape=(n, n, n), extent=((a, b), (a, b), (a, b)))


def test_wavefunction_init_and_basic_properties():
    psi = pt.Wavefunction("exp(-x**2)", g1(128, -1, 1), normalize_const=1.0)
    np.testing.assert_almost_equal(np.sum(psi.state_occupation()), 1.0)

    # multiple internal states, variable usage
    psi2 = pt.Wavefunction(
        ["exp(-x**2)", "b*cos(x)"],
        g1(64, -1, 1),
        variables={"b": 1.0},
        normalize_const=2.0,
    )
    np.testing.assert_almost_equal(np.sum(psi2.state_occupation()), 2.0)

    assert len(psi.r) == 1
    assert len(psi.k) == 1
    assert psi.amp.shape == (128, 1)


def test_wavefunction_dimensional_cases_and_quick_propagation():
    # 1D
    psi_1d_cases = [
        "exp(-x**2)",
        ["exp(-x**2)"],
        ["exp(-x**2)", "0"],
    ]
    V_diag_1d = [
        pt.DiagonalPotential("0"),
        pt.DiagonalPotential("sin(x)"),
        pt.DiagonalPotential("2*x**2"),
    ]

    for initial in psi_1d_cases:
        wf = pt.Wavefunction(initial, g1(32), normalize_const=1.0)
        wf.freely_propagate(steps=1, dt=1.0)

        if wf.num_int_dim == 1:
            for V in V_diag_1d:
                wf.propagate(potential=V, steps=1, dt=1.0)
        else:
            V = pt.HermitianPotential.from_lower_triangular(["0", "sin(x)", "0"])
            wf.propagate(potential=V, steps=1, dt=1.0)

    # 2D
    psi_2d = pt.Wavefunction("exp(-x**2-y**2)", g2(24, 24), normalize_const=1.0)
    psi_2d.propagate(
        potential=pt.DiagonalPotential("x**2 + y**2"),
        steps=1,
        dt=0.1,
    )
    assert psi_2d.amp.shape == (24, 24, 1)

    # 3D + two states, non-diagonal potential
    psi_3d = pt.Wavefunction(["exp(-x**2-y**2-z**2)", "0"], g3(12), normalize_const=1.0)
    psi_3d.propagate(
        potential=pt.HermitianPotential.from_lower_triangular(["0", "x", "sin(y)"]),
        steps=1,
        dt=0.05,
    )


def test_singleton_axes_are_treated_as_inactive_dimensions():
    psi = pt.Wavefunction(
        "exp(-y**2)",
        pt.Grid(shape=(1, 64), extent=((0.0, 0.0), (-2.0, 2.0))),
        normalize_const=1.0,
    )

    assert psi.num_ext_dim == 1
    assert len(psi.r) == 1
    assert psi.r[0].shape == (64,)

    exp_all = psi.exp_pos()
    np.testing.assert_allclose(exp_all, np.array([psi.exp_pos(axis=0)]))

    # Ensure FFT plans are created on active axes only.
    psi.freely_propagate(steps=1, dt=0.01)


def test_state_occupation_conservation_quick():
    wf = pt.Wavefunction(["exp(-x**2)", "0"], g1(64, -3, 3), normalize_const=1.0)
    V = pt.HermitianPotential.from_lower_triangular(["0", "sin(t)", "0"])
    wf.propagate(V, steps=3, dt=0.01)
    np.testing.assert_almost_equal(np.sum(wf.state_occupation()), 1.0, decimal=5)
