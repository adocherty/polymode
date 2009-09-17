from pylab import *
from numpy import *
from Polymode import *

## Parameters:
wl = 1.55
Nx = (200,21)

## Materials
air = Material.Air()
silica = Material.Silica()

r1 = 1
r2 = 2
Daz = 108*pi/180
ring = Waveguide.Annulus(air, r=(r1,r2), phi=(-Daz/2,Daz/2))
wg = Waveguide.Waveguide(material=silica, symmetry=3)
wg.add_shape( ring )

s = NLSolver.DefaultSolver(wg, Nx)
modes = s(wl, 0, 1.4, number=2)
modes += s(wl, 1, 1.4, number=2)

modes.sort(reverse=True)

Plotter.plot_modes_in_grid(modes, 'sz,pe')

show()
