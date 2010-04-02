#!/usr/bin/env python
from os.path import join

#Use setuptools for egg installs, if possible
import setuptools
from numpy.distutils.core import setup, Command

from Polymode import __version__

package_name = 'Polymode'
package_version = __version__
package_description ="A package for the modal analysis of microstructured optical fibers"

class generate_api_docs(Command):
    """Generate the api documentation using epydoc
    """
    description = "generate the api documentation"
    user_options = []

    target_dir = "../documentation/api"

    def initialize_options(self):
        self.all = None

    def finalize_options(self):
        pass

    def run(self):
        import os
        os.system("epydoc --no-frames -o %s Polymode" % self.target_dir)

def configuration(parent_package='',top_path=None):
    from numpy.distutils.misc_util import Configuration

    config = Configuration(None, parent_package, top_path)
    config.set_options(
        ignore_setup_xxx_py=True,
        assume_default_configuration=True,
        delegate_options_to_subpackages=True,
        quiet=False,
    )

    #The module
    config.add_subpackage(package_name)

    #Other packages used
    config.add_subpackage('Nurbs', subpackage_path='other/Nurbs')

    return config

def setup_package():

    setup(
        name = package_name,
        version = package_version,
        description = package_description,
        maintainer = "Andrew Docherty",
        url='http://polymode.googlecode.com',
        license='GPL3',
        configuration = configuration,
#        install_requires = ['numpy >= 1.0.1', 'scipy>=0.5.2', 'matplotlib>=0.92',],
        zip_safe = True,
        cmdclass = {'doc' : generate_api_docs}
        )

    return

if __name__ == '__main__':
    setup_package()

