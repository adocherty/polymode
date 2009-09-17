from __future__ import division

import sys,os
if os.path.pardir not in sys.path:
	sys.path = [os.path.pardir] + sys.path

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname).1s: %(message)s')

from PolyMode.mathlink import timer
timer.reset()

from pylab import *
from numpy import *
from scipy import linalg

from PolyMode import *
from PolyMode import KrylovSolver
reload(KrylovSolver)

#Solver.__polymode_debug__ = True

air = Material.Air()
mat = Material.Fixed(1.45)

wg = Waveguide.Waveguide(rmin=4.25, rmax=9.5, material=mat, symmetry=6)
c = Waveguide.Circle(air, radius=2.0)
wg.shapes = Waveguide.create_hexagonal_tiling(c, layers=1, D=6.75, symmetry=6)

Nx = 200,11					#number of radial & azimuthal points
wl = 1.0							#Wavelength
k0 = 2*pi/wl
m0 = 1
Nmodes = 1

all_convergence = []
all_modes = []

#ksolver = KrylovSolver.NLArnoldiSolver(wg, Nx)
ksolver = KrylovSolver.NLBiArnoldiSolver(wg, Nx)

ksolver(wl, m0, 1.448, number=4)
#ksolver(wl, m0, nefflist=[1.448,1.448])

timer.report()

