from os.path import join
import setuptools

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

package_name = "Nurbs"
version = "0.1.0"

def configuration(parent_package='', top_path=None):
    config = Configuration(parent_package, parent_package, top_path)
    config.add_subpackage('Nurbs')
    return config

if __name__ == '__main__':
    setup(
        name = package_name,
        version = version,
        configuration = configuration,
        zip_safe = True,
        )

