from numpy import *
from Polymode import *

# Materials
core = Material.SiO2GeO2(0.2)
cladding = Material.Silica()

# Create waveguide
wg = Waveguide.Waveguide(material=core, symmetry=1)

#Create the core
s1 = Waveguide.Annulus(cladding, r=(2.0,3.0))
wg.add_shape(s1)

# Waveguide Parameters:
Nx = 100,1				#number of radial & azimuthal points
wl = 1.0					#Wavelength

# Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, number=2)
modes += solver(wl, 1, number=1)

#Plot modes:
Plotter.figure()
Plotter.plot_modes_in_grid(modes, rmax=5)
Plotter.show()

#Save modes
save_data(modes, 'ex2_modes.dat')

