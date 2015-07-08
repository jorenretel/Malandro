#! usr/bin/python

'''This file can be used to only compile c-code and skip the translation
   step from cython to c. The requirement is that this has been done
   before. Typical use:

   python setupC.py build_ext --inplace

   For more installation details, read the README.

'''

from distutils.core import setup
from distutils.extension import Extension
from version import __version__

import numpy
import os.path

current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, 'malandro/backend')
random_dir = os.path.join(current_dir, 'randomGenerator')

randomExtension = Extension("randomGenerator.cyrandom",
                            ["randomGenerator/cyrandom.c"],
                            include_dirs=[numpy.get_include(), src_dir, random_dir])

malandroExtension = Extension("malandro.backend.malandro",
                              ["malandro/backend/malandro.c"],
                              include_dirs=[numpy.get_include(), src_dir, random_dir])

setup(name='Malandro',
      version=__version__,
      description='CCPNMR plug-in for semi-automated sequential assignment.',
      author='Joren Retel',
	  ext_modules=[randomExtension, malandroExtension],
	  packages=['randomGenerator','malandro.gui', 'malandro.backend'])

