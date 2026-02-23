GPU Acceleration
================

pyTALISES supports GPU acceleration via `CuPy <https://cupy.dev/>`_, providing
significant speedups for large-grid simulations.

Installation
------------

Install with GPU support using the optional ``gpu`` extra:

.. code-block:: bash

   pip install pytalises[gpu]

This installs CuPy with CUDA 12.x support. For other CUDA versions, install
CuPy separately:

.. code-block:: bash

   pip install pytalises
   pip install cupy-cuda11x  # For CUDA 11.x

Requirements:

- NVIDIA GPU with CUDA support
- CUDA toolkit (version matching your CuPy install)
- At least one visible CUDA device

Backend Selection
-----------------

By default, pyTALISES auto-detects the best available backend:

.. code-block:: python

   import pytalises as pt

   # Auto-detection: uses CuPy if available, otherwise NumPy
   psi.propagate(potential=V, steps=1000, dt=1e-6)

   # Check what's available
   from pytalises.backends import available_backends, has_cupy
   print(available_backends())  # ('cupy', 'numpy') or ('numpy',)
   print(has_cupy())            # True if CuPy + GPU available

Explicit backend selection via :class:`~pytalises.PropagationOptions`:

.. code-block:: python

   from pytalises import PropagationOptions

   # Force NumPy (useful for debugging or comparison)
   options = PropagationOptions(backend="numpy")
   psi.propagate(potential=V, steps=1000, dt=1e-6, options=options)

   # Force CuPy (raises error if unavailable)
   options = PropagationOptions(backend="cupy")
   psi.propagate(potential=V, steps=1000, dt=1e-6, options=options)

Global default backend:

.. code-block:: python

   from pytalises.backends import set_default_backend

   set_default_backend("cupy")  # All subsequent propagations use CuPy

Performance Characteristics
---------------------------

GPU acceleration is most beneficial for large grids. Typical crossover points
on modern hardware (RTX A4500):

==================  ================  ==============
Grid size           Free propagation  With potential
==================  ================  ==============
n < 2048            NumPy faster      ~Equal
n = 4096            ~Equal            CuPy ~1.2x
n = 8192            CuPy ~1.1x        CuPy ~1.3x
n = 32768           CuPy ~6x          CuPy ~3.5x
n = 131072          CuPy ~12x         CuPy ~14x
==================  ================  ==============

.. note::

   Exact crossover points depend on your specific GPU, CPU, and simulation
   type. The potential workload benefits more at large grids due to the
   analytic 2×2 coupled propagator optimization.

Precision Control
-----------------

Control numerical precision via ``dtype``:

.. code-block:: python

   # Double precision (default) - best accuracy
   options = PropagationOptions(dtype="complex128")

   # Single precision - faster, uses less memory
   options = PropagationOptions(dtype="complex64")

Trade-offs:

- **complex128** (default): Full double precision. Parity errors ~1e-14.
- **complex64**: Single precision. ~1.3-2x faster, but parity errors ~1e-5.

Use ``complex64`` when:

- Speed is critical and moderate precision is acceptable
- GPU memory is constrained
- Running many simulations for parameter sweeps

Advanced Options
----------------

The :class:`~pytalises.PropagationOptions` class provides fine-grained control:

.. code-block:: python

   options = PropagationOptions(
       backend="auto",           # "auto", "numpy", or "cupy"
       dtype="complex128",       # "complex128" or "complex64"
       threads=4,                # CPU threads (NumPy backend)
       coupled_2x2_mode="auto",  # "auto" or "eigh"
       coupled_2x2_kernel="vectorized",  # "vectorized" or "fused"
       profile_stages=False,     # Enable timing profiler
   )

**coupled_2x2_mode**:

- ``auto``: Use analytic closed-form matrix exponential for 2×2 Hermitian
  potentials (much faster, especially on GPU)
- ``eigh``: Force eigendecomposition (slower but handles edge cases)

**coupled_2x2_kernel**:

- ``vectorized``: Standard vectorized operations
- ``fused``: Fused kernel with reduced memory traffic (experimental)

**profile_stages**: When ``True``, timing breakdowns are collected during
propagation. Access via ``psi.stage_timings`` after propagation.

Troubleshooting
---------------

**CuPy not detected:**

.. code-block:: python

   from pytalises.backends import has_cupy, available_backends
   print(has_cupy())           # False
   print(available_backends()) # ('numpy',)

Check:

1. CuPy is installed: ``pip show cupy-cuda12x``
2. CUDA is visible: ``nvidia-smi``
3. GPU is accessible from Python:

   .. code-block:: python

      import cupy as cp
      print(cp.cuda.runtime.getDeviceCount())

**Out of GPU memory:**

- Reduce grid size
- Use ``dtype="complex64"``
- Free other GPU processes

**Results differ between backends:**

Small numerical differences (~1e-14 for complex128) are expected due to
different floating-point ordering. Larger differences may indicate:

- Single precision (complex64) accumulating error over many steps
- Edge cases in potential evaluation

For validation, compare both backends on a short run:

.. code-block:: python

   import numpy as np

   psi_np = psi.copy()
   psi_cp = psi.copy()

   psi_np.propagate(potential=V, steps=100, dt=1e-6,
                    options=PropagationOptions(backend="numpy"))
   psi_cp.propagate(potential=V, steps=100, dt=1e-6,
                    options=PropagationOptions(backend="cupy"))

   diff = np.max(np.abs(psi_np.psi - psi_cp.psi.get()))
   print(f"Max difference: {diff}")  # Should be ~1e-14 for complex128
