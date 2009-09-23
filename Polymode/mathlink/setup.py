#!/usr/bin/env python
from os.path import join

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

#Misc compilation flags
_ublock_compile_with_openmp = False
_ublock_compile = False

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
    
    #Compile the amos library
    config.add_library('famos',
                    sources=[join('amos','*f'),
                    join('amos','mach','*f')]
                    )
    
    #Add the g2c or gfortran libraries to the source, is there a better way?
    def generate_fortran_library(ext, build_dir):
        if config.have_f77c():
            ext.libraries.append('g2c')
        elif config.have_f77c():
            ext.libraries.append('gfortran')
        print "extend libraries:", ext.libraries
        return None

    config.add_extension('bessel_ratios',
                    sources=['bessel_ratios.cpp', generate_fortran_library],
                    depends=['bessel_function_ratios.h'],
                    libraries=['famos']
                    )

    #C++ tridiagonal block solver

    #Get blas/lapack resources from numpy
    lapack_opt = config.get_info('lapack_opt')

    #Get boost SOURCE  resources from numpy
    boost_src = config.get_info('boost_python')

    if not lapack_opt:
        from numpy.distutils.system_info import NotFoundError
        raise NotFoundError,'no lapack/blas resources found'

    ublock_sources = ['ublas_block_lu/ublocklu.cpp','ublas_block_lu/numimport.cpp']
    ublock_depends = ['ublas_block_lu/ublocklu.hpp','ublas_block_lu/numimport.hpp']
    ublock_compile_flags=[]
    ublock_link_flags=[]
    ublock_library_dirs=[]
    ublock_include_dirs=['ublas_block_lu']
    ublock_libraries=[]

    #If boost_src is found, use that, otherwise use default directories
    if boost_src:
        ublock_include_dirs += boost_src['include_dirs']
        ublock_libraries += boost_src['libraries']
    else:
        boost_include_dir = "/usr/include"
        boost_library_dir = "/usr/lib"
        boost_lib_name="boost_python-mt"
        #ublock_library_dirs=[boost_library_dir]
        #ublock_include_dirs=['ublas_block_lu', boost_include_dir]
        #ublock_libraries=['atlas','lapack_atlas',boost_lib_name]
      
    if _ublock_compile_with_openmp:
        ublock_compile_flags += ['-fopenmp', '-ftree-vectorize']
        ublock_link_flags += ['-lgomp']

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
