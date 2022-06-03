#! /usr/bin/env python

# System imports
from distutils.core import *
from distutils      import sysconfig

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
subprocess.run(command.split())

# eec_back extension module
_eec_back = Extension("_eec_back",
                   ["eec_back_wrap.cc", "eec_back.cc"],
                   include_dirs = [numpy_include],
                   libraries=['stdc++', 'm'],
                   extra_compile_args=['-std=c++14']
                   )

# eec_back setup
setup(  name        = "energy-energy correlators",
        description = "tbd",
        author      = "Simon Rothman",
        version     = "0.0",
        ext_modules = [_eec_back]
        )
