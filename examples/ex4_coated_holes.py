from numpy import *
from Polymode import *

## Parameters:
symmetry = 6
wl = 1.45
Nx = 200,21

## Materials
air = Material.Air()
m145 = Material.Fixed(1.45)
m17 = Material.Fixed(1.7)

wg = Waveguide.Waveguide(material=m145, symmetry=symmetry)

#Create a cicular hole filled with air
hole = Waveguide.Circle(air, center=(6.75,0), radius=1.5)

#Create a coating of material m17 from hole with a coating thickness of 0.2um
chole = Waveguide.create_coated_shape(hole, m17, 0.2)

#Add hole and coating to waveguide
wg.add_shape(chole)

#Find mode near 1.45
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 1, 1.45, number=1)

Plotter.figure()
modes[0].plot()
wg.plot(fill=0)
Plotter.show()
