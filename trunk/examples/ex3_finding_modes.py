from pylab import *
from numpy import *
from Polymode import *

# Materials
air = Material.Air()
polymer = Material.Polymer()

# Create waveguide
wg = Waveguide.Waveguide(material=polymer, symmetry=6)

#Construct some shapes to put in the waveguide
wg.add_shape( Waveguide.Circle(air, center=(3,0), radius=1.25) )

# Waveguide Parameters:
Nx = 100,20                                     #number of radial & azimuthal points
wl = 1.0                                                        #Wavelength

# Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)

#Solve for modes in different symmetry classes.
#Notice the "modes +=" adds to the modes
modes = solver(wl, 0, 1.45, number=2)
modes += solver(wl, 1, 1.45, number=2)

# Save the modes and the waveguide
save_data((modes, wg), "ex3_modes.dat")

#Create a new plot window
figure()

Plotter.plot_modes_in_grid(modes)

#This command to show the plot window
#is needed if running outside of ipython
show()
