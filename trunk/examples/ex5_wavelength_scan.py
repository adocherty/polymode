from pylab import *
from numpy import *
from Polymode import *

## Parameters:
sym = 6
Nx = (300,21)

## Materials
air = Material.Air()
pmma = Material.PMMA()

wg = Waveguide.Waveguide(material=pmma, symmetry=sym)

#Create holes in square lattice
layers = 3
D = 6.0
d = 4.0

s = Waveguide.Circle(air, radius=d/2)
shapes = Waveguide.create_hexagonal_tiling(s, layers=layers, D=D, symmetry=sym)
wg.add_shapes(shapes)

#Choose solver
solver = NLSolver.DefaultSolver(wg, Nx)

#Create WL range solver
wlsolver = Solver.WavelengthTrackingSolver(solver)

wlrange=[1.0,1.5]
modes = wlsolver(wlrange, 0, nefflist=[1.3])

Plotter.plot_mode_properties(modes, 'neff', 'wl')

