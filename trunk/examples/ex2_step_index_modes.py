from numpy import *
from Polymode import *

# Materials
core = Material.SiO2GeO2(0.2)
cladding = Material.Silica()

# Create waveguide
wg = Waveguide.Waveguide(material=cladding, symmetry=1)

#Create the core
s1 = Waveguide.Circle(core, center=(0,0), radius=2.0)
wg.add_shape(s1)

# Waveguide Parameters:
Nx = 100,1        #number of radial & azimuthal points
wl = 1.0          #Wavelength

# Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, number=2)
modes += solver(wl, 1, number=1)

for mode in modes:
    mode.info()
