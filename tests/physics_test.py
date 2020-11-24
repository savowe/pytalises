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
        f_Rabi = case[0]
        f_Delta = case[1]
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
