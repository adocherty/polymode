from pylab import *
from numpy import *
from Polymode import *

## Waveguide Parameters:
Nx=200,30
wavelength = 0.764
m0 = 1

## Materials
air = Material.Air()
polymer = Material.Polymer()

## Create waveguide
wg = Waveguide.Waveguide(material=polymer, symmetry=2)

#Construct some shapes to put in the waveguide
radx = 2.1/2; rady = 0.54*radx
Dx = 3.82; yscaling = 0.61

#Construct the waveguide
e = Waveguide.Ellipse(air, axes=(radx,rady))
shapes = Waveguide.create_hexagonal_tiling(e, layers=2, D=Dx, symmetry=2)

#Stretch hexagonal packing
Waveguide.transform_shape_centers(shapes, scale=(1,yscaling))
wg.add_shapes(shapes)

#Create the solvers
allmodes = []
wls = arange(0.1,0.2,0.005)*Dx

birefringence = []
neffapprox = polymer.index(wls[0])-1e-3
for wl in wls:
	solver = NLSolver.DefaultSolver(wg, Nx)
	modes = solver(wl, m0, neffapprox, number=2)

	bifi = 0
	if len(modes)>1:
		neffapprox = modes[0].neff.real
		allmodes += modes

		#Find x polarized mode:
		pa = modes[0].polarization_angle()

		if abs(pa)<pi/4:	#x-polarized
			bifi = modes[1].beta-modes[0].beta
		else:	#y-polarized
			bifi = modes[0].beta-modes[1].beta
		neffapprox = modes[0].neff
	
	birefringence.append(bifi)

Plotter.plot(wls, birefringence, 'b-')
Plotter.xlabel('Wavelength, um')
Plotter.ylabel('Birefringence, ')

