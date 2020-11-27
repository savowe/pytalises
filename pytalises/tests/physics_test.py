import pytalises as pt
import numpy as np
import itertools


def test_rabi_oscillations():
    """
    Tests if the popoulation in a driven two-level system
    undergoes Rabi oscillations. It is compared with analytical
    results.
    """
    f_Rabi_selection = [0.5, 1, 10]
    f_Delta_selection = [0, 0.5, 1]

    cases = itertools.product(f_Rabi_selection, f_Delta_selection)

    for case in cases:
        f_Rabi, f_Delta = case
        f_g_Rabi = np.sqrt(f_Rabi ** 2 + f_Delta ** 2)  # generalized Rabi frequency
        time = 1 / f_g_Rabi

        psi = pt.Wavefunction(["exp(-x**2)", "0"], (128,), (-1, 1), normalize_const=1.0)
        U = pt.Propagator(
            psi,
            ["0", "Omega/2", "Delta"],
            variables={"Omega": 2 * np.pi * f_Rabi, "Delta": 2 * np.pi * f_Delta},
        )

        n_steps = 10
        for _ in range(10):
            np.testing.assert_almost_equal(
                U.psi.state_occupation(),
                np.array(
                    [
                        1
                        - f_Rabi ** 2
                        / f_g_Rabi ** 2
                        * np.sin(2 * np.pi * f_g_Rabi * psi.t / 2) ** 2,
                        f_Rabi ** 2
                        / f_g_Rabi ** 2
                        * np.sin(2 * np.pi * f_g_Rabi * psi.t / 2) ** 2,
                    ]
                ),
            )
            U.potential_prop(delta_t=time / n_steps)
            U.kinetic_prop(delta_t=time / n_steps)


def test_free_propagation_of_gaussian_wave_packet():
    """
    Tests if a gaussian wave packet disperses and moves
    as predicted by the Schrödinger equation.
    """
    sigma0_select = [0.5, 1, 5]
    x0_select = [-3.432, 0, 10]
    k0_select = [-5, 0, 6.2341]
    hbar = 1.054571817e-34
    m_select = [hbar, hbar * 2.2]

    cases = itertools.product(sigma0_select, x0_select, k0_select, m_select)

    for case in cases:
        print(case)
        sigma0, x0, k0, m = case
        psi = pt.Wavefunction(
            ["exp(-((x-x0)/(2*sigma0))**2)*exp(1j*k0*x)"],
            (32768,),
            (-200, 200),
            normalize_const=1.0,
            variables={"sigma0": sigma0, "x0": x0, "k0": k0},
            m=m,
        )

        for _ in range(10):
            np.testing.assert_almost_equal(
                psi.var_pos(),
                sigma0 ** 2 + (hbar ** 2 * psi.t ** 2) / (4 * m ** 2 * sigma0 ** 2),
                decimal=2,
            )
            np.testing.assert_almost_equal(
                psi.exp_pos(), x0 + hbar / m * k0 * psi.t, decimal=2
            )
            psi.freely_propagate(num_time_steps=1, delta_t=1)


def test_two_dimensional_harmonic_oscillator():
    """
    Tests if the center of mass of a gaussian wave packet in
    an harmonic oscillator moves as expected.
    """

    x0_select = [0, 0.5]
    y0_select = [0.5, 1.5]
    omega_x_select = [2 * np.pi * 0.8]
    omega_y_select = [0, 2 * np.pi * 1]

    cases = itertools.product(x0_select, y0_select, omega_x_select, omega_y_select)

    for case in cases:
        x0, y0, omega_x, omega_y = case
        psi = pt.Wavefunction(
            ["exp(-((x-x0)/(2*1))**2)*exp(-((y-y0)/(2*1))**2)"],
            (128, 128),
            [(-7.5, 7.5), (-7.5, 7.5)],
            normalize_const=1.0,
            variables={"x0": x0, "y0": y0},
        )
        time = 2 * np.pi / omega_x
        n_steps = 10
        for _ in range(n_steps):
            psi.propagate(
                "1/2*omega_x**2*x**2 + 1/2*omega_y**2*y**2",
                variables={"omega_x": omega_x, "omega_y": omega_y},
                diag=True,
                num_time_steps=10,
                delta_t=1 / 10 * time / n_steps,
            )
            np.testing.assert_almost_equal(
                psi.exp_pos(),
                np.array([x0 * np.cos(omega_x * psi.t), y0 * np.cos(omega_y * psi.t)]),
                decimal=1,
            )
