from pylab import *
from numpy import *
from Polymode import *

import pickle

## Waveguide Parameters:
Nx=100,11               #number of radial & azimuthal points

## Materials
air = Material.Air()
m2 = Material.Fixed(2.0)
polymer = Material.Polymer()

## Create waveguide
wg = Waveguide.Waveguide(material=polymer, symmetry=6)

#Construct some shapes to put in the waveguide
wg.add_shape(Waveguide.Circle(air, center=(3,0), radius=1.25))

#Create the solver
neffstart = 1.4
wl = 1.45
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, neffstart, number=2)

#set compress_mode flag to size of mode information to store
solver = NLSolver.DefaultSolver(wg, Nx, compress_to_size = (10,5))
modes_comp = solver(wl, 0, neffstart, number=2)

#set compress_mode flag to size of mode information to store
solver = NLSolver.DefaultSolver(wg, Nx, store = False)
modes_novec = solver(wl, 0, neffstart, number=2)

print "Size of uncompressed modes: %dkb" % (len(pickle.dumps(modes))/1024)
print "Size of compressed modes: %dkb" % (len(pickle.dumps(modes_comp))/1024)
print "Size of modes without vectors: %dkb" % (len(pickle.dumps(modes_novec))/1024)

clf()
subplot(121)
modes[0].plot()

subplot(122)
modes_comp[0].plot()
