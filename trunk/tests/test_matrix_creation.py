from __future__ import division

import sys,os
if os.path.pardir not in sys.path:
	sys.path = [os.path.pardir] + sys.path

import logging
logging.basicConfig(level=logging.DEBUG, format='%(levelname).1s: %(message)s')

from pylab import *
from numpy import *

from PolyMode import *
from PolyMode.mathlink import *
from PolyMode.difflounge import boundary,finitedifference
from PolyMode.mathlink import blockarray

Nx = (7,1)
bcwidth = 3
bandwidth = 3
xbc = 0.5
m0 = 1
wl = 1.0

mat = Material.Fixed(1.45)
wg = Waveguide.Waveguide(rmin=2.5, rmax=9.5, material=mat, symmetry=1)

coord = wg.get_coord(Nx, border=0)

#Boundary conditions are Neumann for m+m0==0 and Dirichlet otherwise
bcl = boundary.BC_Mixed(coord.rv, bandwidth=bcwidth, xbc=xbc, dtype=complex_)
#bcl = boundary.BC_OnePoint(coord.rv, dtype=complex_)
bcr = boundary.BC_Mixed(coord.rv, bandwidth=bcwidth, xbc=xbc, dtype=complex_)
#bcr = boundary.BC_OnePoint(coord.rv, dtype=complex_)

#Setup default finite differences
fd_diff = finitedifference.DifferenceMatrix(2, bandwidth=bandwidth, \
			X=coord.rv, bcr=bcr, bcl=bcl, dtype=complex_)
fd_jac = finitedifference.JacobianMatrix(2, bandwidth=bandwidth, \
			X=coord.rv, bcr=bcr, bcl=bcl, dtype=complex_)

eq = Equation.VectorWaveEquation(fd_diff)
jac = Equation.VectorWaveJacobian(fd_jac)

eq.setup(Nx,wg,m0,wl)
jac.setup(Nx,wg,m0,wl)

ev = (1.4275+1e-6j)**2*eq.k0**2

eq.set_lambda(ev)

#Construct analytic solution
n1 = wg.exterior_index(wl)
g = sqrt(n1**2*eq.k0**2 - ev+0j)
rv = coord.calc_rv(1)

Nr, Naz = Nx
bw = eq.diff.bandwidth
pmax = eq.pmax
blockshape = (pmax*Naz,)*2
A = blockarray.BlockArray((Nr,bw), blockshape=blockshape, dtype=complex_)
M = A.blockview()

#Extents over which the boundary condition changes the lines
bcstart,bcend = eq.diff.bc_extents()
rows = range(Nr)
vec_row = zeros((1,Naz,pmax), dtype=complex_)
for ii in rows:
	blockstart = max(bw//2 - ii,0)
	for kk in range(pmax*Naz):
		vec_row.flat[kk] = 1
		y = eq.construct(vec_row, ii)
		M[ii,blockstart:blockstart+y.shape[0],:,kk] = y.reshape((y.shape[0], 2*Naz))
		vec_row.flat[kk] = 0

Nmax = pmax*Naz*Nr
A_direct = zeros((Nmax,)*2, dtype=complex_)
vec_row = zeros((Nr,Naz,pmax), dtype=complex_)
for ii in range(Nmax):
	vec_row.flat[ii] = 1
	y = eq.rmatvec(vec_row)
	A_direct[:,ii] = y.flat
	vec_row.flat[ii] = 0

A_direct = A_direct.T.conj()

print A_direct.real
print A.toarray().real
print "Error", abs(A.toarray()-A_direct).max()

