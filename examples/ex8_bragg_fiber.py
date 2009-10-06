##
## Bragg fiber
##
from numpy import *
from pylab import *

from Polymode import *
from Polymode import LayeredSolver

ncore = 1.0
n1=1.4
n2=1.6

rcore = 5.0
Nlayer = 20
m = 0

neff_opt = 0.998
wl_opt = 1.0

d1 = wl_opt/(4*sqrt(n1**2 - neff_opt**2))
d2 = wl_opt/(4*sqrt(n2**2 - neff_opt**2))
dlayer = d1+d2

wlrange=[0.4,1.4]

print "Bragg layer widths d1=%.4g, d2=%.4g, rcore = %.5g " % (d1,d2, rcore)

mg = Material.Fixed(ncore)
m1 = Material.Fixed(n1)
m2 = Material.Fixed(n2)

#create structure object with external substrate of m1
wg = Waveguide.Waveguide(material=m1, symmetry=1)

#Core
wg.add_shape(Waveguide.Circle(mg, center=(0,0), radius=rcore))

#Cladding
for ii in range(Nlayer):
    ri = rcore+ii*dlayer
    a1 = Waveguide.Annulus(m1, r=(ri,ri+d1))
    a2 = Waveguide.Annulus(m2, r=(ri+d2,ri+dlayer))
    wg.add_shapes(a1,a2)

#Use the LayeredSolver
solver = LayeredSolver.LayeredSolver(wg)

#Do an overall scan of the condition number of the eigenvalue
#problem
wlrange=[0.8,1.2]
wlcscan = Solver.WavelengthConditionScan(solver, Nscan=(50,50))
wlcscan(wlrange, m, neffrange=[ncore-0.05, ncore])
wlcscan.plot(style={'cmap':cm.gray_r})

#Track the fundamental mode 
wlrange=[0.95,1.15]
wltrack = Solver.WavelengthTrack(solver)
wltrack.ga_target=1e-4
wltrack(wlrange, m, nefflist=[0.992], number=1)
wltrack.plot()

