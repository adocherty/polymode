#!/usr/bin/env python
from os.path import join

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

try:
    import Cython.Compiler.Main
    have_cython = True
except ImportError:
    have_cython = False

#Fix build_src to work with Cython
from numpy.distutils.command import build_src

if have_cython and not build_src.have_pyrex:
	print "Numpy doesn't support Cython!"

def configuration(parent_package='', top_path=None):

    config = Configuration('mathlink', parent_package, top_path)
    config.add_library('famos',sources=[join('amos','*f'), join('amos','mach','*f')])
    #config.add_extension('amos', sources=['amos.cpp'], depends=['amos.h'], libraries=['famos', 'g2c'])
    config.add_extension('bessel_ratios', sources=['bessel_ratios.cpp'], depends=['bessel_function_ratios.h'], libraries=['famos', 'gfortran'])

	#C++ tridiagonal block solver
    #config.add_extension('ublocklu', sources=['ublas_block_lu/ublocklu.cpp','ublas_block_lu/numimport.cpp'], libraries=[])

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
