from setuptools import setup


setup(name='pytalises',
      version='0.1.1',
      description="""
      TALISES (This Ain't a LInear Schrödinger Equation Solver) is an easy-to-use Python implementation
      of the Split-Step Fourier Method, for numeric calculation of a wave function's time-propagation
      under the Schrödinger equation.""",
      url='http://github.com/savowe/pytalises',
      author='Sascha Vowe',
      author_email='sascha.vowe@posteo.de',
      license='GPL v3.0',
      packages=['pytalises'],
      install_requires=[
          'numpy',
          'scipy',
          'numba',
          'pyfftw',
          'numexpr',
      ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "Development Status :: 3 - Alpha",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics"
    ],
    python_requires='>=3.6',
    zip_safe=False)

