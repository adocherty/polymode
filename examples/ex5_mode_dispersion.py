# _*_ coding=utf-8 _*_
from numpy import *
from Polymode import *

#Basic materials
air = Material.Air()
m145 = Material.Fixed(1.45)
pmma = Material.Polymer()

Nx = (200,41)
wl = 1.0
dwl = 1e-3

#The waveguide consists of 6 holes spaces 6.75um apart in a hexagon
wg = Waveguide.Waveguide(material=pmma, symmetry=6)
wg.add_shape( Waveguide.Circle(air, center=(6.75,0), radius=2.5) )

#Calculate modes at three different wavelengths
solver = NLSolver.DefaultSolver(wg, Nx)
modes = []
for wlx in [wl-dwl, wl, wl+dwl]:
    modes += solver(wlx, 1, number=1)

#Estimate dβ/dλ with a finite difference
dbetadwl = (modes[2].beta - modes[0].beta)/(2*dwl)

m = modes[1]

#Calculate on full symmetry
ng = -dbetadwl*wl**2/(2*pi)
ngm = m.group_index(wg)

print "\nCalculations of group index:\n"
print "ng from wavelength:",ng
print "ng from integral:", ngm
