from numpy import *
from Polymode import *

## Waveguide Parameters:
Nx=400,60               #number of radial & azimuthal points
wl = 1.4
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
fn_index = lambda d: 1-(d/2)**4.25

#The low index cladding - put it below other shapes
cladding = Waveguide.Circle(mlow, center=(0,0), radius=rcladding, zorder=-1)

#The high index core
core = Waveguide.Circle(silica, center=(0,0), radius=radius1)

#The high index inclusions - put the high index dopant above other shapes
rings = []
for i in range(Nrings):
    rings += [ Waveguide.Circle(silica, center=((i+1)*D,0), radius=radius1) ]
    inclusion = Waveguide.Circle(mhigh, center=((i+1)*D,0), radius=radius2, zorder=1)
    inclusion.set_index_function(fn_index, background=silica)
    rings += [ inclusion ]

wg.add_shapes(core, cladding, rings)

bifi_last = None
Nxs = [(400,20), (400,40), (400,80), (400,160)]
Nxs = [(200,20)]
for Nx in Nxs:
    #Create the solver
    solver = NLSolver.DefaultSolver(wg, Nx)

    #Find the two fundamental core modes (we hope)
    modes = solver(wl, m0, nefflist=[1.4429]*2)

    bifi = modes[0].beta-modes[1].beta
    print "Birefringence = %.5g at Nx=%s" % (bifi.real, Nx)

    if bifi_last:
        print "Change:", abs(bifi-bifi_last)
    bifi_last = bifi

Plotter.plot_modes_in_grid(modes, 'sz', cartesian=1)
Plotter.plot_modes_in_grid(modes, 'vectore', wg=wg)

