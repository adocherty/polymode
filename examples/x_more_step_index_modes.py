from __future__ import division
from numpy import *

from Polymode import *
from Polymode import LayeredSolver

# Materials
core = Material.Fixed(1.48)
cladding = Material.Silica()

# Create waveguide
wg = Waveguide.Waveguide(material=cladding, symmetry=1)

#Create the core
s1 = Waveguide.Circle(core, center=(0,0), radius=3.0)
wg.add_shape(s1)

# Waveguide Parameters:
Nx = 500,1                      #number of radial & azimuthal points
wl = 1.0                                        #Wavelength

# Create the solvers
wlrange=[1.0,3.0]

solver = NLSolver.DefaultSolver(wg, Nx)
wlsolver = Solver.WavelengthTrackingSolver(solver)
modes = wlsolver(wlrange, 1,  neffrange=[cladding.index(wl), core.index(wl)], number=3)

solver = LayeredSolver.DefaultSolver(wg)
wlsolver = Solver.WavelengthTrackingSolver(solver)
modesl = wlsolver(wlrange, 1,  neffrange=[cladding.index(wl), core.index(wl)], number=3)

Plotter.plot_mode_properties(modes, 'neff', 'wl', 'bx')
Plotter.plot_mode_properties(modesl, 'neff', 'wl', 'r.')
