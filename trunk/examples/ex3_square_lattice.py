from pylab import *
from numpy import *
from Polymode import *

## Parameters:
sym = 4
Nx = (200,40)
wl = 1.0

## Materials
air = Material.Air()
si = Material.Silica()

wg = Waveguide.Waveguide(rmax=6.4, material=si, symmetry=sym)

#Create holes in square lattice
layers = 2
D = 2.5
d = 0.5*D

c = Waveguide.Circle(air, radius=d/2)
shapes = Waveguide.create_square_tiling(c, layers=layers, Dx=D, Dy=D, symmetry=sym)
wg.add_shapes(shapes)

#Find modes
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 0, neffrange=[1.4,1.45], number=5)
modes += solver(wl, 1, neffrange=[1.4,1.45], number=5)
modes += solver(wl, 2, neffrange=[1.4,1.45], number=5)

Plotter.figure()
Plotter.plot_modes_in_grid(modes, wg=wg)
Plotter.show()
