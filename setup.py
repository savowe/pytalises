from setuptools import setup

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(
    name="pytalises",
    version="0.2.1",
    description="""
      TALISES (This Ain't a LInear Schrödinger Equation Solver) is an easy-to-use Python implementation
      of the Split-Step Fourier Method, for numeric calculation of a wave function's time-propagation
      under the Schrödinger equation.""",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="http://github.com/savowe/pytalises",
    author="Sascha Vowe",
    author_email="sascha.vowe@posteo.de",
    license="GPL v3.0",
    packages=["pytalises"],
    install_requires=[
        "numpy>=1.19",
        "scipy>=1.5",
        "numba>=0.50",
        "pyfftw>=0.12",
        "numexpr>=2.7",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    python_requires=">=3.6",
    zip_safe=False,
)
