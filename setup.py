#! /usr/bin/env python

# System imports
from setuptools import setup, Extension

# Third-party modules - we depend on numpy for everything
import numpy

# Obtain the numpy include directory.  This logic works across numpy versions.
try:
    numpy_include = numpy.get_include()
except AttributeError:
    numpy_include = numpy.get_numpy_include()

import subprocess

command = 'swig -python -c++ -fastproxy -keyword -py3 -o eec_back_wrap.cc eec_back.i'
print(command)
subprocess.run(command.split(), cwd='eec/backend/')

# eec_back extension module
_eec_back = Extension("_eec_back",
                   ["eec/backend/eec_back_wrap.cc", "eec/backend/eec_back.cc"],
                   include_dirs = [numpy_include],
                   libraries=['stdc++', 'm'],
                   extra_compile_args=['-std=c++14']
                   )

# eec_back setup
setup(  name        = "energy-energy correlators",
        description = "tbd",
        author      = "Simon Rothman",
        version     = "0.1.0",
        ext_modules = [_eec_back],
        packages = ['eec']
        )
