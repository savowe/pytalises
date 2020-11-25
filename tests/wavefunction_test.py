import pytalises as pt
import numpy as np
import itertools


def test_Wavenfunction_init():
    # 1D Arguments not as list
    psi = pt.Wavefunction("exp(-x**2)", (128,), (-1, 1), normalize_const=1.0)
    np.testing.assert_almost_equal(psi.state_occupation(), 1.0)
    # 1D Arguments as list
    pt.Wavefunction(["exp(-x**2)"], (128,), [(-1, 1)])
    pt.Wavefunction(["exp(-x**2)", "cos(x)"], (128,), [(-1, 1)])
    pt.Wavefunction(["exp(-x**2)", "b"], (128,), [(-1, 1)], variables={"b": 1})
    # 2D Arguments not as list
    pt.Wavefunction("exp(-x**2-y**2)", (128, 128), [(-1, 1), (-1, 1)])
    # 2D Arguments as list
    pt.Wavefunction(["exp(-x**2-y**2)"], (128, 128), [(-1, 1), (-1, 1)])
    psi = pt.Wavefunction(
        ["exp(-x**2-y**2)", "c*x**2+y**2"],
        (128, 128),
        [(-1, 1), (-1, 1)],
        variables={"c": 1},
        normalize_const=666,
    )
    np.testing.assert_almost_equal(np.sum(psi.state_occupation()), 666)
    # 3D Arguments not as list
    pt.Wavefunction(
        "exp(-x**2-y**2-z**2/d)+d",
        (128, 128, 128),
        [(-1, 1), (-1, 1), (-1, 1)],
        variables={"d": 1},
    )
    # 3D Arguments as list
    pt.Wavefunction(
        ["exp(-x**2-y**2-z**2/d)+d"],
        (128, 128, 128),
        [(-1, 1), (-1, 1), (-1, 1)],
        variables={"d": 1},
    )

    # 1D + one internal state
    psi = [
        "exp(-x**2)",
        ["exp(-x**2)"],
    ]
    number_of_grid_points = [
        16,
        (16,),
        (
            16,
            1,
        ),
        (16, 1, 1),
    ]
    spatial_ext = [
        (-2, 2),
        [
            (-2, 2),
        ],
    ]
    variables = [
        {"a0": 0, "x0": 1},
    ]
    V = ["0", ["0"], "sin(x)", "2*x**2"]
    all_valid_cases = itertools.product(
        psi, number_of_grid_points, spatial_ext, variables, V
    )
    for case in all_valid_cases:
        psi = pt.Wavefunction(
            psi=case[0],
            number_of_grid_points=case[1],
            spatial_ext=case[2],
            variables=case[3],
        )
        psi.freely_propagate(num_time_steps=1, delta_t=1)
        psi.propagate(potential=case[4], num_time_steps=1, delta_t=1)

    # 2D + one internal state
    psi = [
        "exp(-x**2)*y",
        ["exp(-x/a0**2)*exp(-y**2)"],
    ]
    number_of_grid_points = [
        (
            16,
            16,
        ),
        (16, 16, 1),
    ]
    spatial_ext = [
        [(-2, 2), (-5, 5)],
    ]
    variables = [
        {"a0": 0, "x0": 1},
    ]
    V = ["0", ["0"], "sin(x*y)", "2*x**2*y"]
    all_valid_cases = itertools.product(
        psi, number_of_grid_points, spatial_ext, variables, V
    )
    for case in all_valid_cases:
        psi = pt.Wavefunction(
            psi=case[0],
            number_of_grid_points=case[1],
            spatial_ext=case[2],
            variables=case[3],
        )
        psi.freely_propagate(num_time_steps=1, delta_t=1)
        psi.propagate(potential=case[4], num_time_steps=1, delta_t=1)

    # 3D + two internal state
    psi = [
        ["exp(-x/a0**2)*exp(-y**2)*exp(-z**2)", "0"],
    ]
    number_of_grid_points = [
        (16, 16, 16),
    ]
    spatial_ext = [
        [(-2, 2), (-5, 5), (-2, 2)],
    ]
    variables = [
        {"a0": 0, "x0": 1},
    ]
    V = [
        ["0", "sin(x*y)", "z"],
        ["x", "sin(x*y)"],
    ]
    all_valid_cases = itertools.product(
        psi, number_of_grid_points, spatial_ext, variables, V
    )
    for case in all_valid_cases:
        if len(case[4]) == 2:
            diag = True
        else:
            diag = False
        psi = pt.Wavefunction(
            psi=case[0],
            number_of_grid_points=case[1],
            spatial_ext=case[2],
            variables=case[3],
        )
        psi.freely_propagate(num_time_steps=1, delta_t=1)
        psi.propagate(potential=case[4], num_time_steps=1, delta_t=1, diag=diag)
