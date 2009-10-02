from pylab import *
from numpy import *
from Polymode import *

## Waveguide Parameters:
Nx=100,20
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

solver = NLSolver.DefaultSolver(wg, Nx)

wlrange = array([0.1,0.2])*Dx
neffapprox = polymer.index(wlrange[0])-1e-4

#Track fundamental modes over wavelength range
wltrack = Solver.WavelengthTrack(solver, dont_lose_modes=True)
wltrack.ga_target=1e-3
modes = wltrack(wlrange, m0, neffapprox, number=2)

#Sort modes by polarization and calculate birefringence
plot_wls=[]
plot_bifi=[]

wls = array([m.wl for m in modes])
for wl in wls:
    wlinx = nonzero(wls==wl)
    mx = array(modes)[wlinx]
    if len(mx)==2:
        #Find x polarized mode:
        pa = mx[0].polarization_angle()

        if abs(pa)<pi/4: #x-polarized
            bifi = mx[1].beta - mx[0].beta
        else:                #y-polarized
            bifi = mx[0].beta - mx[1].beta

        plot_bifi.append(bifi)
        plot_wls.append(wl)

Plotter.plot(plot_wls, plot_bifi, 'b-')
Plotter.xlabel('Wavelength, um')
Plotter.ylabel('Birefringence, ')

