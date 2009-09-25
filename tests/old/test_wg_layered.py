from __future__ import division

import sys,os
if os.path.pardir not in sys.path:
	sys.path = [os.path.pardir] + sys.path

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname).1s: %(message)s')

from pylab import *
from numpy import *

from PolyMode import *
from PolyMode import LayeredSolver

from PolyMode.mathlink import timer

ng,n1,n2=1.0,1.46,1.2

wl = 1.1
dlayer = 0.5
alayer = (1+sqrt(n1**2-1)/sqrt(n2**2-1))**(-1)
rcore = 5.0
Nlayer = 10
m = 0

mg = Material.Fixed(ng)
mg.color="white"
m1 = Material.Fixed(n1)
m1.color="lightblue"
m2 = Material.Fixed(n2)
m2.color="darkblue"

#Construct the layers
ns = [ng] + [n1,n2]*Nlayer
gammas = [0]+[0,0]*Nlayer

#Construct layers
wg = Waveguide.Waveguide(material=m2)
wg.add_shape( Waveguide.Circle(mg, radius=rcore) )
for ii in range(Nlayer):
	r1 = rcore + dlayer*ii
	wg.add_shape( Waveguide.Annulus(m1, r=(r1,r1+alayer*dlayer)) )
	wg.add_shape( Waveguide.Annulus(m2, r=(r1+alayer*dlayer,r1+dlayer)) )

#Find all modes
#lsolver = LayeredSolver.LayeredSolver(wg)
#lsolver(m, wl, neffrange=[0.6,1.0], nefflist=[0.999], number=inf)

#Nx = (1000,1)
#solver = Solver.DefaultSolver(Nx, wg)
#solver(m, wl, 0.999, number=1)


