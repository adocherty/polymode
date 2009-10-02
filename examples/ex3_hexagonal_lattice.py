from pylab import *
from numpy import *
from Polymode import *

## Parameters:
sym = 6
Nx = (200,30)
wl = 1.0

## Materials
air = Material.Air()
pmma = Material.PMMA()

wg = Waveguide.Waveguide(material=pmma, symmetry=sym)

#Create holes in square lattice
layers = 2
D = 6.0
d = 4.0

s = Waveguide.Circle(air, radius=d/2)
shapes = Waveguide.create_hexagonal_tiling(s, layers=layers, D=D, symmetry=sym)
wg.add_shapes(shapes)

#Find modes
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, nefflist=[1.48, 1.475,1.46])

Plotter.figure()
Plotter.plot_modes_in_grid(modes, wg=wg)
Plotter.show()
