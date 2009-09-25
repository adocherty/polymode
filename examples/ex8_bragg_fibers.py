##
## Bragg fibers
##
from numpy import *
from pylab import *

from Polymode import *
from Polymode import LayeredSolver
reload(LayeredSolver)

ncore = 1.44
n1=1.46
n2=1.44
neff_opt = 1.436
wl_opt = 1.0

d1 = wl_opt/(4*sqrt(n1**2 - neff_opt**2))
d2 = wl_opt/(4*sqrt(n2**2 - neff_opt**2))
dlayer = d1+d2

#Core size
rcore = 0.5*(3*d1 + 4*d2)
Nlayer = 20
m = 0
wlrange=[0.4,1.4]

print "Bragg layer widths d1=%.4g, d2=%.4g, rcore = %.5g " % (d1,d2, rcore)

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
    a1 = Waveguide.Annulus(m1, r=(ri,ri+d1))
    a2 = Waveguide.Annulus(m2, r=(ri+d2,ri+dlayer))
    wg.add_shapes(a1,a2)

#Create WL range solver using standard FD method
Nx=(1000,1)

#The default FD solver
#solver = NLSolver.DefaultSolver(wg, Nx)

solver = LayeredSolver.LayeredSolver(wg)
solver.setup(remove_boundary_factors = True, debug_plot=True)

#Layered solver using the transfer matrix method
solver = LayeredSolver.LayeredSolverCauchy(wg)
solver.setup(remove_boundary_factors = True, Nscan=10, debug_plot=True)

modes = []
wls = arange(wlrange[0], wlrange[1], (wlrange[1]-wlrange[0])/50)
wls=[1.0]
for wl in wls:
    modes += solver(wl, m, neffrange=[ncore-0.005,ncore],  number=5)

#The wavelength solver takes another solver as the only argument
#wlsolver = Solver.WavelengthTrackingSolver(solver)

#Guess accuracy needs to be more accurate for tracking with closely spaced modes
#wlsolver.ga_target = 1e-5
#modes = wlsolver(wlrange, m, neffrange=[ncore-0.1,ncore],  number=2)

if len(modes)>1:
    subplot(211)
    Plotter.plot_mode_properties(modes, 'neff', 'wl', 'g.')
    subplot(212)
    Plotter.plot_mode_properties(modes, 'loss', 'wl', 'g.')
    Plotter.semilogy()
    #draw()
