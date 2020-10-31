[![PyPI](https://img.shields.io/pypi/v/pytalises?color=blue)](https://pypi.org/project/pytalises/)
[![Conda Version](https://img.shields.io/conda/vn/conda-forge/pytalises?color=blue&label=conda-forge)](https://anaconda.org/conda-forge/pytalises)
[![Build Status](https://travis-ci.com/savowe/pytalises.svg?token=nZF2LbDmAxqpxqs5m7HE&branch=master)](https://travis-ci.com/savowe/pytalises)
[![Documentation Status](https://readthedocs.org/projects/pytalises/badge/?version=latest)](https://pytalises.readthedocs.io/en/latest/?badge=latest)
# pyTALISES

**TALISES** (This Ain't a LInear Schrödinger Equation Solver) is an easy-to-use Python implementation of the Split-Step Fourier Method, for numeric calculation of a wave function's time-propagation under the Schrödinger equation.

### Documentation
Read the [documentation](https://pytalises.readthedocs.io/en/latest/) to learn more about pytalises' capabilities.

### Features
- Calculation of a wavefunction's time propagation under a (non)linear Schrödinger equation: ![](https://latex.codecogs.com/png.latex?%5Cdpi%7B120%7D%20i%5Chbar%20%5Cfrac%7B%5Cpartial%7D%7B%5Cpartial%20t%7D%20%5CPsi%20%28%5Cvec%7Br%7D%2C%20t%29%20%3D%20%5CBig%5BV%28%5CPsi%2C%5Cvec%7Br%7D%2C%20t%29%20&plus;%20%5Cfrac%7B%5Chbar%5E2%7D%7B2m%7D%5Cnabla%5E2%20%5CBig%5D%20%5CPsi%20%28%5Cvec%7Br%7D%2C%20t%29)
- the wave-function ![](https://latex.codecogs.com/gif.latex?%5Cdpi%7B100%7D%20%5CPsi) may include an arbitrary number of internal and external degrees of freedom
- simple implementation of Hamiltonians
- speed of the [FFTW](https://pypi.org/project/pyFFTW/), [BLAS](https://www.netlib.org/blas/) and [numexpr](https://numexpr.readthedocs.io/en/latest/) libaries with multithreading
- crucial functions are just-in-time compiled with [numba](https://numba.readthedocs.io/en/stable/)

Installing pytalises
====================
**We recommend installing pytalises via conda**

#### Using conda

Installing `pytalises` from the `conda-forge` channel can be achieved by adding `conda-forge` to your channels with:

```
conda config --add channels conda-forge
```

Once the `conda-forge` channel has been enabled, `pytalises` can be installed with:

```
conda install pytalises
```


#### Using pip

pytalises is available on the Python Package Index and can be installed via

```
pip install pytalises
```

it has dependencies via `scipy` and `numba` on BLAS and LAPACK libraries that are not always found on windows systems. For linux they can usually be located.
