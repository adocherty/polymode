from numpy import *
from Polymode import *

# Materials
air = Material.Air()
mco = Material.Fixed(1.5)
mcl = Material.Fixed(1.45)

#Create waveguide
wg = Waveguide.Waveguide(material=mcl, symmetry=2)

#The core
s1 = Waveguide.Circle(mco, center=(0,0), radius=5)
wg.add_shape(s1)

#Waveguide Parameters:
Nx=400,1                                                #number of radial & azimuthal points
wavelength = 1.55               #Wavelength
m0 = 1                                                  #Symmetry class

# Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wavelength, m0, number=2)

for m in modes:
    m.normalize(by='ext')

    print
    print "Mode effective index:", m.neff
    print "Loss: %.3g dB" % m.loss
    print "Symmetry class:", m.m0
    print "Propagation constant:", m.beta
    print "Wavelength:", m.wl

    print "Group index:", m.group_index(wg)
    print "Propagation contant from integral:", m.integral_propagation_lossless(wg)
    print "Proportion of power in core:", real(m.mode_power(r=wg.core_size)/m.mode_power())
