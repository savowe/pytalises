V1 to V2 Migration Guide
========================

pyTALISES v2 introduces a cleaner public API around explicit grids,
structured potentials, and consolidated runtime options.

Quick mapping
-------------

- ``number_of_grid_points`` + ``spatial_ext`` -> :class:`pytalises.Grid`
- ``diag=True/False`` + string lists ->
  :class:`pytalises.DiagonalPotential` / :class:`pytalises.HermitianPotential`
- ``num_time_steps`` -> ``steps``
- ``delta_t`` -> ``dt``
- ``num_of_threads`` / ``FFTWflags`` -> :class:`pytalises.PropagationOptions`
- optional GPU execution -> ``PropagationOptions(backend="cupy")`` with
  ``Wavefunction(..., backend="cupy")``

Backend internals
-----------------

The backend/engine internals are private implementation details in v2.
Only the public simulation API is considered stable:
:class:`pytalises.Grid`, :class:`pytalises.Wavefunction`,
:func:`pytalises.propagate`, :func:`pytalises.freely_propagate`, and
structured potentials/options.

Example (before)
----------------

.. code-block:: python

   import pytalises as pt

   psi = pt.Wavefunction(
       ["exp(-x**2)", "0"],
       (256,),
       (-4, 4),
       variables={"x0": 0},
   )

   psi.propagate(
       potential=["V0*x**2", "Omega*cos(t)", "V1*x**2"],
       num_time_steps=1000,
       delta_t=1e-6,
       diag=False,
       variables={"V0": 1.0, "V1": 1.0, "Omega": 2.0},
   )

Example (after)
---------------

.. code-block:: python

   import pytalises as pt

   grid = pt.Grid(shape=(256,), extent=((-4, 4),))
   psi = pt.Wavefunction(
       initial=["exp(-x**2)", "0"],
       grid=grid,
       variables={"x0": 0},
   )

   V = pt.HermitianPotential.from_lower_triangular(
       ["V0*x**2", "Omega*cos(t)", "V1*x**2"]
   )

   psi.propagate(
       potential=V,
       steps=1000,
       dt=1e-6,
       variables={"V0": 1.0, "V1": 1.0, "Omega": 2.0},
       options=pt.PropagationOptions(threads=4),
   )

Legacy adapter
--------------

A temporary compatibility adapter is available as ``pytalises.legacy``.

.. code-block:: python

   import pytalises as pt

   psi = pt.legacy.Wavefunction("exp(-x**2)", (256,), (-4, 4))
   psi.freely_propagate(num_time_steps=100, delta_t=1e-4)

Using ``pytalises.legacy`` emits ``DeprecationWarning`` and will be removed in a future release.
