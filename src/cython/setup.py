#! usr/bin/python

from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
import numpy
import os.path

current_dir = os.path.dirname(os.path.abspath(__file__))
src_dir = os.path.join(current_dir, 'randomGenerator')


randomExtension = Extension("randomGenerator.cyrandom",
                            ["randomGenerator/cyrandom.pyx"],
                            include_dirs = [numpy.get_include(), current_dir, src_dir])

malandroExtension = Extension("malandro",
                            ["malandro.pyx"],
                            include_dirs = [numpy.get_include(), current_dir, src_dir])

randomTestExtension = Extension("randomTest",
                            ["randomTest.pyx"],
                            include_dirs = [numpy.get_include(), current_dir, src_dir])

setup(
    cmdclass = {'build_ext': build_ext},
    ext_modules = cythonize([randomExtension, malandroExtension, randomTestExtension]),
    #ext_modules = cythonize([randomExtension, randomTestExtension]),
    packages=["malandroPackage", "malandroPackage.randomGenerator"]
    
    #[Extension("malandro", ["malandro.pyx"]), Extension("cyrandom", ["cyrandom.pyx"])],
    #include_dirs = [numpy.get_include(),],
)
