#!/usr/bin/env python
from __future__ import division

#--------------------------------------------------------------------------
#Add ABCSolver path as parent directory -- not needed when installed
import sys, os
sys.path.append(os.path.pardir)
#--------------------------------------------------------------------------

from numpy import *
from ABCSolver import *

#Basic materials
air = Material.Air()
m145 = Material.Fixed(1.45)

## Parameters:
symmetry = 6
Nx = (100,11)
wl = 1.0

wg = Waveguide.Waveguide(material=m145, symmetry=symmetry)
wg.add_shape( Waveguide.Circle(air, center=(6.75,0), radius=2.5) )

#Calculate modes at two different wavelengths
solver = NLSolver.DefaultSolver(Nx, wg)
modes = []
for wlx in [wl-0.01, wl+0.01]:
	modes += solver(1, 2*pi/wlx, [1.45, 1.4], number=1)
	

