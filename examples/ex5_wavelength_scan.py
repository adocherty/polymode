from pylab import *
from numpy import *
from Polymode import *

## Parameters:
sym = 2
Nx = (200,20)

## Materials
air = Material.Air()
pmma = Material.PMMA()

wg = Waveguide.Waveguide(material=pmma, symmetry=sym)

layers = 1
D = 6.0
d = 4.0
wl=1.0

s = Waveguide.Rectangle(air, axes=(d,d))
shapes = Waveguide.create_hexagonal_tiling(s, layers=layers, D=D, symmetry=sym)
wg.add_shapes(shapes)

#Choose solver
solver = NLSolver.DefaultSolver(wg, Nx)

#modes=solver(wl, 1, neffrange=[1.47, 1.478], number=4)

#Create WL range solver
wlsolver = Solver.WavelengthTrack(solver)

wlrange = [1.0,1.5]
modes = wlsolver(wlrange, 1, nefflist=[1.48028,1.48028], number=1)

subplot(121)
Plotter.plot_mode_properties(modes, 'neff', 'wl')
subplot(122)
Plotter.plot_mode_properties(modes, 'loss', 'wl')

