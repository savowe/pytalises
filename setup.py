from setuptools import setup


setup(name='pytalises',
      version='0.1.1',
      description='',
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
      zip_safe=False)

