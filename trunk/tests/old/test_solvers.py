#test_waveguide.py

#--------------------------------------------------------------------------
#Remove main ABCSolver pathname & add current path to use the local version
import sys, os, logging
sys.path.append("../../")
logging.basicConfig(level=logging.INFO, format='%(levelname).1s: %(message)s')
#--------------------------------------------------------------------------
from pylab import *

from ABCSolver import Waveguide, Solver, Material
reload(Solver)

air = Material.Air()
pm = Material.Polymer()

Nx = 400,41
wg = Waveguide.Waveguide(9.4, material=pm, symmetry=6, oversample=10)
wg.add_shape( Waveguide.Circle(air, center=(6.75,0), radius=2.5) )

k0=2*pi/1.45
m0=1

#solver_classes = [CenterSolver, BracketSolver]
bs = Solver.BracketSolve(Nx, wg)
cs = Solver.CenterSolve(Nx, wg)

#Check direct call
print "Checking %s\n" % bs.__class__.__name__
bs(m0, k0, [1.48,1.4], number=20, restart=5)

#Check deffered run
#svr.initialize(m0, k0, [1.45,1.1], totalnumber=2)
#svr.calculate(1)
#svr.calculate(1)

#Check direct call
print "Checking %s\n" % cs.__class__.__name__
cs(m0, k0, [1.48,1.4], intervals=4, number=20, restart=5)

#close = Solver.ClosestSolve(Nx, wg)
#close(m0, k0, [1.45, 1.4, 1.3])

print "bs modes:"
for m in sorted(bs.modes):
	print m
	
print "cs modes:"
for m in sorted(cs.modes):
	print m
	
