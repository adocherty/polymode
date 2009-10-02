from numpy import *
from Polymode import *

# Core radius and number of points
Nx = 500,1
k = 2*pi/1.0

## Materials
air = Material.Air()
silica = Material.Silica()
dsilica = Material.SiO2GeO2(0.2)        #20% doped silica

## Create waveguide
wg = Waveguide.Waveguide(rmax=8, material=silica, symmetry=2)

## Parabolic index profile based on distance from center of object
fn_parabolic = lambda d: 1-(d/5)**2

## Doped core
s1 = Waveguide.Circle(dsilica, center=(0,0), radius=5)
s1.set_index_function(fn_parabolic)

wg.add_shape(s1)

Plotter.figure()
wg.plot_bitmap(Nx, k0=k)
Plotter.show()
