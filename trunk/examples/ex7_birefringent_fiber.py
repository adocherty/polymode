from numpy import *
from Polymode import *

## Waveguide Parameters:
Nx=200,41               #number of radial & azimuthal points
wl = 1.18
m0 = 1

## Materials
silica = Material.Silica()
high = Material.SiO2GeO2(0.2)   #Silica with 20% Germania
clad = Material.SiO2Fl(0.011)           #Silica with 1.1% Flourine

## Create waveguide
wg = Waveguide.Waveguide(material=clad, symmetry=2)

#Construct some shapes to put in the waveguide
radius1 = 3.5
radius2 = 3
D = 6.8
Nrings = 2

#Refractive index function
fn_parabolic = lambda d: 1-(d/3)**2

#Construct the waveguide
core = Waveguide.Circle(silica, center=(0,0), radius=radius1)
rings = []
for i in range(Nrings):
    rings += [ Waveguide.Circle(silica, center=((i+1)*D,0), radius=radius1) ]
    inclusion = Waveguide.Circle(high, center=((i+1)*D,0), radius=radius2, zorder=1)
    inclusion.set_index_function(fn_parabolic, background=silica)
    rings += [ inclusion ]

wg.add_shapes(core, rings)

#Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, m0, number=2)

Plotter.plot_modes_in_grid(modes, 'sz', cartesian=1)
Plotter.plot_modes_in_grid(modes, 'vectore', wg=wg)
