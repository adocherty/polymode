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

Nx = (50,1)
bcwidth = 3
bandwidth = 3
xbc = 0.5
m0 = 1
wl = 1.0

mat = Material.Fixed(1.45)
wg = Waveguide.Waveguide(rmin=2.0, rmax=9.5, material=mat, symmetry=1)

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

ev = (1.4475+1e-9j)**2*eq.k0**2

#Construct analytic solution
n1 = wg.exterior_index(wl)
g = sqrt(n1**2*eq.k0**2 - ev+0j)
rv = coord.calc_rv(1)

h = zeros(Nx+(2,), dtype=complex_)
for ii in range(len(rv)//2):
	h[ii] = jv(eq.msm0+eq.pv, g*rv[ii])

jfac = 1 #jv(eq.msm0+eq.pv, g*rv[ii])/hankel1(eq.msm0+eq.pv, g*rv[ii])
for ii in range(len(rv)//2,len(rv)):
	h[ii] = hankel1(eq.msm0+eq.pv, g*rv[ii]) * jfac

#Setup the equation
eq.set_lambda(ev)

#Apply equation
err = abs(eq.matvec(h) - ev*h)
err = err/err.max()

plot(rv,err[:,0,0],'-', rv[-1:], err[-1:,0,0],'r.', rv[:1], err[:1,0,0],'rx')
plot(rv,err[:,0,1],'--', rv[-1:],err[-1:,0,1],'b.', rv[:1],err[:1,0,1],'bx')
semilogy()
draw()

rdtn = rv[-1]
dtn = g*hankel1p(eq.msm0+eq.pv, g*rdtn)/hankel1(eq.msm0+eq.pv, g*rdtn)
print "Right DtN match: ", abs(dtn-eq.diff.bcright.dtn).max()

rdtn = rv[0]
dtn_left = g*jvp(eq.msm0+eq.pv, g*rdtn)/jv(eq.msm0+eq.pv, g*rdtn)
print "Left DtN match: ", abs(dtn_left-eq.diff.bcleft.dtn).max()
print

print "Right DtN error"
print abs( h[-1] * dtn*2*coord.dr + h[-2] - hankel1(eq.msm0+eq.pv, g*coord.rv[-1]) * jfac)
print err[-1]
print

print "Left DtN error"
print abs( -h[0] * dtn_left*2*coord.dr + h[1] - jv(eq.msm0+eq.pv, g*coord.rv[0]) * jfac)
print err[0]
print

#h[:]=0;h[0]=1

#Check Jacobian
print "Jacobian:"
jac.set_lambda(ev)
dh =jac.matvec(h)

dev=1e-5
eq.set_lambda(ev+dev)
dh_num = eq.matvec(h)
eq.set_lambda(ev-dev)
dh_num -= eq.matvec(h)
dh_num /= 2*dev

print "Right", abs(dh[-3:]-dh_num[-3:]).max()
print "Left:", abs(dh[:3]-dh_num[:3]).max()
print

print "Adjoint Jacobian:"
jac.set_lambda(ev)
dh =jac.rmatvec(h)

dev=1e-5
eq.set_lambda(ev+dev)
dh_num = eq.rmatvec(h)
eq.set_lambda(ev-dev)
dh_num -= eq.rmatvec(h)
dh_num /= 2*dev

print "Right", abs(dh[-3:]-dh_num[-3:]).max()
print "Left:", abs(dh[:3]-dh_num[:3]).max()

