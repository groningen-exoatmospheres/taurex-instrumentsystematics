import setuptools
from setuptools import find_packages
from setuptools import setup
#from numpy.distutils.core import Extension
#from numpy.distutils import log
import re, os

packages = find_packages(exclude=('tests', 'doc'))
provides = ['taurex_instrumentsystematics', ]

requires = []

version = "1.0.0" 

install_requires = ['taurex']

entry_points = {'taurex.plugins': 'taurex_instrumentsystematics = taurex_instrumentsystematics'}

setup(name='taurex_instrumentsystematics',
      author="Quentin Changeat",
      author_email="q.changeat@rug.nl",
      license="BSD",
      description='Add multiple spectra and Instrument Systematics fitting to TauREx 3.1',
      packages=packages,
      version=version,
      url="https://github.com/groningen-exoatmospheres/taurex-instrumentsystematics",
      entry_points=entry_points,
      provides=provides,
      requires=requires,
      install_requires=install_requires)
