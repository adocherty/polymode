from numpy import *
from Polymode import *

# Materials
air = Material.Air()
m2 = Material.Fixed(2.0)
polymer = Material.Polymer()

# Create waveguide
wg = Waveguide.Waveguide(material=polymer, symmetry=6)

#Construct some shapes to put in the waveguide
s1 = Waveguide.Ellipse(air, center=(3,pi/6), axes=(1.0,1.4), rot=pi/6)
s2 = Waveguide.RegularPolygon(air, center=(6,0), Nsides=6, length=2.0, rot=pi/6)
s3 = Waveguide.RegularPolygon(m2, center=(6,0), Nsides=6, length=1.5, rot=pi/6, zorder=1)

wg.add_shapes(s1, s2, s3)

wl = 1.45                                               #Wavelength
Nx = 100,21                             #number of radial & azimuthal points

# Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)

#Solve for modes in different symmetry classes.
#here we give a starting neff to search near
modes = solver(wl, 1, 1.449, number=4)

#Find mode with lowest losses:
losses = [m.loss for m in modes]
llmode = modes[argmin(losses)]

print "Found low loss mode:", llmode

Plotter.figure()

#Plot the mode
llmode.plot()

#Plot the waveguide overlayed
wg.plot(fill=0)

Plotter.show()
