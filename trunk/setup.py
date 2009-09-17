#!/usr/bin/env python
from os.path import join
import setuptools

package_name = 'Polymode'
version = '0.1.0'

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=True,
    )
    
    #The module
    config.add_subpackage(package_name)

    #Other packages used
    config.add_subpackage('Nurbs', subpackage_path='other/Nurbs')

    return config

def setup_package():

    from numpy.distutils.core import setup

    setup(
        name = package_name,
        version = version,
        maintainer = "Andrew Docherty",
        maintainer_email = "docherty@gmail.com",
        url='http://polymode.googlecode.com',
        license='GPL3',
        configuration = configuration,
        install_requires = ['numpy >= 1.0.0', 'scipy>=0.7.0', 'matplotlib>0.90',],
        zip_safe = True,
        )

    return

if __name__ == '__main__':
    setup_package()

