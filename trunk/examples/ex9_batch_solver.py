from pylab import *
from numpy import *
from Polymode import *

## Parameters:
symmetry = 4
Nx = (200,21)
wl = 1.0

## Materials
air = Material.Air()
si = Material.Silica()

wg = Waveguide.Waveguide(material=si, symmetry=symmetry)

#Create holes in square lattice
Nrings = 2
D = 4.0
d = 0.5*D
for nr in range(1,Nrings+1):
    for ii in range(Nrings+1):
        wg.add_shape( Waveguide.Circle(air,center=(nr*D,ii*D),radius=d/2,xy=1) )

aff = pi/(2*sqrt(3))*(d/D)**2
print "Created waveguide with %d rings and air fill fraction f=%.2g" % (Nrings, aff)

#Set up batched solver
m0s = [0,1,2]

bsolve = Solver.BatchSolve(NLSolver.DefaultSolver, wg, Nx, compress_to_size=(25,10))
bsolve.initialize(m0s, wl, totalnumber=100)
bsolve.set_neffrange([1.3, 1.45], intervals=10)

modes = Solver.batch_solve(bsolve)

"""
from numpy import arange,sqrt, random, linalg
from multiprocessing import Pool

global counter
counter = 0
def cb(r):
    global counter
    print counter, r
    counter +=1

def det(M):
    return linalg.det(M)

po = Pool()
for i in xrange(1,300):
    j = random.normal(1,1,(100,100))
    po.apply_async(det,(j,),callback=cb)
po.close()
po.join()
print counter
"""

print tick.report()
