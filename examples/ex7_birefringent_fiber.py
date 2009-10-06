from numpy import *
from Polymode import *

## Waveguide Parameters:
Nx=300,61               #number of radial & azimuthal points
wl = 1.1
m0 = 1

## Materials
silica = Material.Silica()
mhigh = Material.SiO2GeO2(0.25)           #Silica with 20% Germania
mlow = Material.SiO2Fl(0.011)           #Silica with 1.1% Flourine

## Create waveguide
wg = Waveguide.Waveguide(material=silica, symmetry=2)

#Construct some shapes to put in the waveguide
D = 8.0
radius1 = D/2
radius2 = D/4
Nrings = 4

rcladding = (2*Nrings+1)*radius1

#Refractive index function
fn_parabolic = lambda d: 1-(d/3)**2

#The low index cladding - put it below other shapes
cladding = Waveguide.Circle(mlow, center=(0,0), radius=rcladding, zorder=-1)

#The high index core
core = Waveguide.Circle(silica, center=(0,0), radius=radius1)

#The high index inclusions
rings = []
for i in range(Nrings):
    rings += [ Waveguide.Circle(silica, center=((i+1)*D,0), radius=radius1) ]
    inclusion = Waveguide.Circle(mhigh, center=((i+1)*D,0), radius=radius2, zorder=1)
    inclusion.set_index_function(fn_parabolic, background=silica)
    rings += [ inclusion ]

wg.add_shapes(core, cladding, rings)

#Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)

#Find the two modes closest to the RI of the core
modes = solver(wl, m0, 1.447,  number=10)

Plotter.plot_modes_in_grid(modes, 'sz', cartesian=1)
Plotter.plot_modes_in_grid(modes, 'vectore', wg=wg)
