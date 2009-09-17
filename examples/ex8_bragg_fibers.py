##
## Bragg fibers
##
from numpy import *
from pylab import *

from Polymode import *
from Polymode import LayeredSolver

ncore = 1.0
n1=1.6
n2=1.4
rcore = 5.0
dlayer = 0.5
Nlayer = 20
wlrange=[0.8,1.4]

alayer = (1+sqrt(n1**2-1)/sqrt(n2**2-1))**(-1)
m = 0

mg = Material.Fixed(ncore)
m1 = Material.Fixed(n1)
m2 = Material.Fixed(n2)

#create structure object with external substrate of m2
wg = Waveguide.Waveguide(material=m2, symmetry=1)

#Core
wg.add_shape(Waveguide.Circle(mg, center=(0,0), radius=rcore))

#Cladding
for ii in range(Nlayer):
	ri = rcore+ii*dlayer
	a1 = Waveguide.Annulus(m1, r=(ri,ri+alayer*dlayer))
	a2 = Waveguide.Annulus(m2, r=(ri+alayer*dlayer,ri+dlayer))
	wg.add_shapes(a1,a2)

#Create WL range solver using standard FD method
Nx=(1000,1)

solver = NLSolver.DefaultSolver(wg, Nx)
wlsolver = Solver.WavelengthTrackingSolver(solver)
modes = wlsolver(wlrange, 0, nefflist=[0.995,0.995])

#Create WL range solver using the transfer matrix method
#solver = LayeredSolver.DefaultSolver(wg)
#wlsolver = Solver.WavelengthTrackingSolver(solver)
#lmodes = wlsolver(wlrange, 0, nefflist=[0.995,0.995])

subplot(211)
Plotter.plot_parameters(modes, 'neff', 'wl', 'g.')
#Plotter.plot_parameters(lmodes, 'neff', 'wl', 'r-')

subplot(212)
Plotter.plot_parameters(modes, 'loss', 'wl', 'g.')
#Plotter.plot_parameters(lmodes, 'loss', 'wl', 'r-')

draw()
