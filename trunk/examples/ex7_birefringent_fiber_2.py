from pylab import *
from numpy import *
from Polymode import *

## Solver parameters
Nx=200,41
m0 = 1

## Waveguide Parameters
radius1 = 3.5
radius2 = 3
D = 6.8
Nrings = 2

## Materials
silica = Material.Silica()
high = Material.SiO2GeO2(0.2)   #Silica with 20% Germania
clad = Material.SiO2Fl(0.011)           #Silica with 1.1% Flourine

## Create waveguide
wg = Waveguide.Waveguide(material=clad, symmetry=2)

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

#Solve at difference wavelengths
wls=arange(1.18,1.6,0.005)

modes=[]
allmodes=[]
birefringence=[]
bifi=0
for wl in wls:
    if len(modes)>1:
        modes = solver(wl, m0, modelist=modes)
    else:
        neffapprox=silica.index(wl)
        modes = solver(wl, m0, neffapprox, number=2)

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

