import pytalises as pt


def test_Wavenfunction():
    # 1D Arguments not as list
    pt.Wavefunction("exp(-x**2)", (128,), (-1, 1))
    # 1D Arguments as list
    pt.Wavefunction(["exp(-x**2)"], (128,), [(-1, 1)])
    pt.Wavefunction(["exp(-x**2)", "cos(x)"], (128,), [(-1, 1)])
    pt.Wavefunction(["exp(-x**2)", "b"], (128,), [(-1, 1)], variables={"b": 1})
    # 2D Arguments not as list
    pt.Wavefunction("exp(-x**2-y**2)", (128, 128), [(-1, 1), (-1, 1)])
    # 2D Arguments as list
    pt.Wavefunction(["exp(-x**2-y**2)"], (128, 128), [(-1, 1), (-1, 1)])
    pt.Wavefunction(
        ["exp(-x**2-y**2)", "c*x**2+y**2"],
        (128, 128),
        [(-1, 1), (-1, 1)],
        variables={"c": 1},
    )
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
