#!/usr/bin/env python
from os.path import join

from numpy.distutils.core import setup
from numpy.distutils.misc_util import Configuration

def configuration(parent_package='', top_path=None):

    config = Configuration('Polymode', parent_package, top_path)
    config.add_subpackage('mathlink')
    config.add_subpackage('difflounge')
    config.add_data_dir('data')
    return config

if __name__ == '__main__':
    setup(**configuration(top_path='').todict())
