#Test the difference on a 2d array
from numpy import *
from scipy import sparse, linalg

from ABCSolver.mathlink import timer

from ABCSolver.difflounge import finitedifference, linesparse, boundary

from ABCSolver.difflounge.finitedifference import *
from ABCSolver.difflounge.boundary import *
from ABCSolver.difflounge.linesparse import *


#Timing tests
tick = timer.timer()

Nx=1000; Ny=1
X=20.; Y=1.
dx=X/Nx; dy=Y/Ny

v = random.random((Nx)).astype(complex_)

#The complete (including boundary points) X,Y space
x = arange(0,Nx+2)*dx; y = arange(0,Ny)*dy

bcleft = BC_Mixed(X=x, bandwidth=3, condition=1, xbc=0, dtype=complex_)
bcright = BC_Mixed(X=x, bandwidth=3, condition=1, xbc=0, dtype=complex_)

#fd = SparseDifferenceMatrix(2, bandwidth=5, X=x, bcl=bcleft, bcr=bcright, dtype=float)
#fdm = fd.generate()

t = timer.timer()
t.clear()

fdas = DifferenceMatrix(2, bandwidth=5, X=x, bcl=bcleft, bcr=bcright, dtype=complex_)
fdas.generate()
D1 = fdas.dmatrix[1].toarray()

t.start('dot')
for ii in range(10):
	yd = dot(D1,v)
t.stop('dot')

t.start('fdm')
for ii in range(10):
	y = fdas.diff1(v)
t.stop('fdm')

print v.shape, y.shape

print "Difference", abs(yd-y).max()
t.report()

