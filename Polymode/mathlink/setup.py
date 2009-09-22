#!/usr/bin/env python
from os.path import join

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

#Misc compilation flags
_ublock_compile_with_openmp = True
_ublock_compile = False

#boost_include_dir = "/opt/boost/include"
#boost_library_dir = "/opt/boost/lib"
#boost_lib_name="boost_python-gcc43-mt"
boost_include_dir = "/usr/include"
boost_library_dir = "/usr/lib"
boost_lib_name="boost_python-mt"

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

    ublock_sources = ['ublas_block_lu/ublocklu.cpp','ublas_block_lu/numimport.cpp']
    ublock_depends = ['ublas_block_lu/ublocklu.hpp','ublas_block_lu/numimport.hpp']
    ublock_library_dirs=[boost_library_dir]
    ublock_include_dirs=['ublas_block_lu', boost_include_dir]
    ublock_libraries=['atlas','lapack_atlas',boost_lib_name]
    ublock_compile_flags=[]
    ublock_link_flags=[]

    if _ublock_compile_with_openmp:
        ublock_compile_flags += ['-fopenmp', '-ftree-vectorize']
        ublock_link_flags += ['-lgomp']

    #C++ tridiagonal block solver
    if _ublock_compile:
      config.add_extension('ublocklu',
                    sources=ublock_sources,
                    depends=ublock_depends,
                    library_dirs=ublock_library_dirs,
                    include_dirs=ublock_include_dirs,
                    libraries=ublock_libraries,
                    extra_compile_args=ublock_compile_flags,
                    extra_link_args=ublock_link_flags,
                    extra_info = lapack_opt
                    )

    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
