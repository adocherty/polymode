import pylab as pl
from numpy import *

from Polymode import *
from Polymode.mathlink import blockarray

air = Material.Air()
pmma = Material.Polymer()
pmma = Material.Fixed(1.49)

## Parameters:
symmetry = 6
m0 = 1
wl = 1
Nx = 300,50

sw = 0.4                #Wall width
s1 = 7.4                #Primary scaling parameter: the length of the inner walls
Nt = 2                  #Number of layers

#Calculate wall length of dodecagons
s2 = s1/(1+2*cos(pi/6))

#Tiling distance
dt = s1*cos(pi/6) + s1/2

#Build waveguide from exact shapes:
wg = Waveguide.Waveguide(material=pmma, symmetry=symmetry)

#Larger central core hole:
p1 = Waveguide.RegularPolygon(air, Nsides=6, length=s1, zorder=2)

#Coat the air hole with the thin wall of pmma width sw
shapes = Waveguide.create_coated_shape(p1, pmma, sw, outer=True, inner=True)

#Cladding as a tiled dodecagon, again coating an air hole with pmma
p2 = Waveguide.RegularPolygon(air, Nsides=12, length=s2)
p2c = Waveguide.create_coated_shape(p2, pmma, sw, outer=True, inner=True)
shapes += Waveguide.create_hexagonal_tiling(p2c, layers=Nt, D=dt, symmetry=symmetry)

#Put a background of air so the interstitial triangles are now air-filled
shapes.append(Waveguide.RegularPolygon(air, Nsides=6, length=2*s2*sqrt(3)*Nt, rot=pi/6, zorder=-2))

wg.add_shapes(shapes)

#Solve for some modes
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, m0, 1.0, number=6)

Plotter.plot_modes_in_grid(modes, cartesian=1)
