from __future__ import division
from pylab import *
from numpy import *

from Polymode import *
from Polymode import Optimize

## Waveguide Parameters:
Nx=200,30
wl = 1.0
m0 = 1

## Materials
air = Material.Air()
polymer = Material.Polymer()

#Number of ellipses:
Ne = 6

def adjust(s, p):
    ## Create waveguide
    shapes = []

    ps = reshape(p, (Ne,3))

    #Construct the waveguide
    for ax,cx,cy in ps:
        shapes.append(Waveguide.Ellipse(air, center=(cx,cy), axes=(ax,1.5*ax), rot=0, xy=1))

    s.wg.shapes = shapes

def evaluate(s):
    #Optimize for the loss of the lowest loss mode
    neffs = [m.neff for m in s.modes]

    return imag(neffs).min()

wg = Waveguide.Waveguide(material=polymer, symmetry=2)
solver = NLSolver.DefaultSolver(wg, Nx)
solver.initialize(wl, m0, neffrange=[1.42, 1.48], number=3)

xs = [2.5,2.5,5,7.5,7.5,0]
ys = [-5,5,0,-4,4,7.5]

pstart = []
for ii in range(Ne):
    pstart += [1.0, xs[ii], ys[ii]]

o=Optimize.Optimize(solver, adjust, evaluate)
out = o(pstart)
