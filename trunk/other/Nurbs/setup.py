from os.path import join

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):

    config = Configuration('Nurbs', parent_package, top_path)
    config.add_extension('_Bas', sources=['_Bas.c'])
    config.add_data_dir('Doc')

    return config

if __name__ == '__main__':
    setup(version='0.1',
        description='Nurbs - Python module to work with NURBS curves and surfaces',
        author='Runar Tenfjord',
        author_email = 'runten@netcom.no',
        url='http://runten.tripod.com/',
        zip_safe=True,
        **configuration(top_path='').todict())
