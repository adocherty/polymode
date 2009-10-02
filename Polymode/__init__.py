# _*_ coding=utf-8 _*_
"""
Polymode
========

The Polymode Modal Solver is a project is a library for the calculation
of mode patterns and parameters in optical fibers, especially microstructured 
optical fiber (MOF) or photonic crystal fibers (PCF) with silica or polymer 
construction.

The focus of Polymode is to find large numbers of leaky modes automatically in
microstructured optical fibers of arbitrary cross section using novel mode
calculation algorithms.

This project intends to create a library to find modes of optical fibers with 
arbitrary cross-section robustly. The code is based on Python and Numpy/Scipy,
with extentions in C++. Matplotlib is required for the plotting functionality.

-----

Copyright © 2009 Andrew Docherty

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""
__version__ = '0.1'

# Other modules that should be loaded for 'from Polymode import *':
__all__ = ['Modes', 'Image', 'Solver', 'NLSolver', 'Waveguide', 'Material', 'Plotter', 'Equation']
__all__ += ['save_data', 'load_data', 'coordinates', 'constants']

#Easy configure logging
import logging
logging.basicConfig(level=logging.INFO, format='%(levelname).1s: %(message)s')

#Bring to top level
from .mathlink import coordinates, constants

#Some helper functions
def save_data(data, filename="data.mds"):
    from cPickle import dump
    f = open(filename, "wb")
    dump(data, f)
    f.close()

def load_data(filename="data.mds"):
    from cPickle import load
    f = open(filename, "rb")
    data = load(f)
    f.close()
    return data

