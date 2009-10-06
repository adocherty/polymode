from pylab import *
from numpy import *
from Polymode import *

## Solver parameters
Nx=300,61
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

#Solve at difference wavelengths
wls=arange(1.12,1.6,0.005)
neffapprox = 1.4467

modes=[]
allmodes=[]
birefringence=[]
for wl in wls:
    bifi = nan

    if len(modes)>1:
        modes = solver(wl, m0, modelist=modes)
    else:
        modes = solver(wl, m0, nefflist=[neffapprox, neffapprox])

    if len(modes)>1:
        #Find x polarized mode:
        pa = modes[0].polarization_angle()

        if abs(pa)<pi/4:        #x-polarized
            bifi = modes[0].beta-modes[1].beta
        else:   #y-polarized
            bifi = modes[1].beta-modes[0].beta
        neffapprox = modes[0].neff

    print "B:",bifi
    birefringence.append(bifi)
    allmodes.extend(modes)

#Compress modes
allmodes = Modes.compress_modes(allmodes, (20,11), wg)

#Save all data
save_data((wls, birefringence, allmodes), "ex7_bifi.dat")

plot(wls, birefringence)
xlabel(r'Wavelength, $\mu$m')
ylabel(r'Birefringence, $\beta_x-\beta_y$')

