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
    from numpy.distutils.system_info import get_info, NotFoundError

    config = Configuration('mathlink', parent_package, top_path)
    
    #Compile the amos library
    config.add_library('famos',
                    sources=[join('amos','*f'),
                    join('amos','mach','*f')]
                    )

    config.add_extension('bessel_ratios',
                    sources=['bessel_ratios.cpp'],
                    depends=['bessel_function_ratios.h'],
                    libraries=['famos', 'gfortran']
                    )

    #Get blas/lapack resources from numpy
    lapack_opt = get_info('lapack_opt')

    if not lapack_opt:
        raise NotFoundError,'no lapack/blas resources found'

    atlas_version = ([v[3:-3] for k,v in lapack_opt.get('define_macros',[]) \
                      if k=='ATLAS_INFO']+[None])[0]

    #C++ tridiagonal block solver
    ublock_sources = ['ublas_block_lu/ublocklu.cpp','ublas_block_lu/numimport.cpp']
    config.add_extension('ublocklu',
                    sources=ublock_sources,
                    include_dirs=['ublas_block_lu'],
                    libraries=['atlas','lapack_atlas','boost_python-mt'],
                    extra_info = lapack_opt
                    )

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
