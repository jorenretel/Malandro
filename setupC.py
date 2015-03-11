#! usr/bin/python

'''This file can be used to only compile c-code and skip the translation
   step from cython to c. The requirement is that this has been done
   before. Typical use:

   python setup.py build_ext --inplace

   For more installation details, read the README.

'''

from distutils.core import setup
from distutils.extension import Extension

import numpy
import os.path

current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, 'src/cython')
random_dir = os.path.join(current_dir, 'src/cython/randomGenerator')

randomExtension = Extension("src.cython.randomGenerator.cyrandom",
                            ["src/cython/randomGenerator/cyrandom.c"],
                            include_dirs=[numpy.get_include(), src_dir, random_dir])

malandroExtension = Extension("src.cython.malandro",
                              ["src/cython/malandro.c"],
                              include_dirs=[numpy.get_include(), src_dir, random_dir])

setup(ext_modules=[randomExtension, malandroExtension])

