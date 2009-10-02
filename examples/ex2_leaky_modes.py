from numpy import *
from Polymode import *
from Polymode import LayeredSolver

# Materials
core = Material.SiO2GeO2(0.2)
cladding = Material.Silica()

# Create waveguide
wg = Waveguide.Waveguide(material=core, symmetry=1)

#Create the core
s1 = Waveguide.Annulus(cladding, r=(2.0,3.0))
wg.add_shape(s1)

# Waveguide Parameters:
Nx = 100,1                              #number of radial & azimuthal points
wl = 1.0                                        #Wavelength

# Create the solver
#solver = LayeredSolver.LayeredSolver(wg)
solver = NLSolver.DefaultSolver(wg, Nx)

#Solve for 2 modes in mode class 0
modes = solver(wl, 0, number=2)

#Solve for 1 mode in mode class 1
modes += solver(wl, 1, number=1)

#Plot modes, this defaults to 1d, as we specify Nphi = 1
#We plot to a radial distance of 5
Plotter.plot_modes_in_grid(modes, rmax=5)
Plotter.show()

#Save modes
save_data(modes, 'ex2_modes.dat')
