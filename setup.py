#! usr/bin/python

'''Setup file. This file can be used to compile the cython source
   to a working python extension. To do this, both cython and
   a c compiler like gcc should be installed. Typical use:

   python setup.py build_ext --inplace

   For more installation details, read the README.

'''

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import os.path

current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, 'malandro/backend')
random_dir = os.path.join(current_dir, 'randomGenerator')

randomExtension = Extension("randomGenerator.cyrandom",
                            ["randomGenerator/cyrandom.pyx"],
                            include_dirs=[numpy.get_include(), src_dir, random_dir])

malandroExtension = Extension("malandro.backend.malandro",
                              ["malandro/backend/malandro.pyx"],
                              include_dirs=[numpy.get_include(), src_dir, random_dir])

setup(cmdclass={'build_ext': build_ext},
      ext_modules=cythonize([randomExtension, malandroExtension]),
      packages=['randomGenerator', 'malandro.gui', 'malandro.backend'])
