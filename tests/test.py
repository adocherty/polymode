#Test the difference on a 2d array
from __future__ import division

from numpy import *
from scipy import sparse, linalg

from ABCSolver.mathlink import timer

from ABCSolver.difflounge.finitedifference import *
from ABCSolver.difflounge.boundary import *
from ABCSolver.difflounge.line_sparse import *

#Timing tests
tick = timer.timer()

Nx=400; Ny=2
X=2.; Y=1.
dx=X/Nx; dy=Y/Ny

v = random.random((Nx,Ny))

#The complete (including boundary points) X,Y space
x = arange(0,Nx+2)*dx; y = arange(0,Ny)*dy

#Construct sparse arrays
Fx = sparse.lil_matrix((Nx, Nx), dtype=complex)

#Construct this order of difference
Nd = 2

#Finite difference bandwidth
bw = 7
offset=bw//2

#boundary condition stencil width
bcbw = 5
bcorder = 1
bcoffset = bcbw//2
xint = x[bcorder:-bcorder]

#Boundary conditions
#bcleft = BC_Dirichlet()
#bcright = BC_Dirichlet()
m1=0.5
m2=-0.5
bcleft = BC_Mixed(X=x, bandwidth=bcbw, condition=m1, xbc=0, dtype=complex)
bcright = BC_Mixed(X=x, bandwidth=bcbw, condition=m2, xbc=1, dtype=complex)
#bcleft = BC_OnePoint(X=x, condition=m1)
#bcright = BC_OnePoint(X=x, condition=m2)

#Extents over which the boundary condition changes the lines
bcstart = offset
bcend = bw-offset-1
print "Boundary condition affects rows: 0->%d and %d->N" % (bcstart, bcend)
tick.start("sparse, main")
#Set main elements
for ii in range(bcstart,Nx-bcend):
	iis = max(0,ii-offset); iie = min(ii+bw-offset,Nx)
	Fx[ii, iis:iie] = fd_calc_weights(xint[iis:iie],xint[ii],Nd)[Nd]
tick.stop("sparse, main")

tick.start("sparse, BC")
#Left boundary conditions
for ii in range(0,bcstart):
	xiis = 0; xiie = ii+bw-offset+bcorder
	
	#Apply boundary conditions to stencil based on x with boundary
	#nodes included, but centered at ii+bcoffset
	w = fd_calc_weights(x[xiis:xiie],xint[ii],Nd)[Nd]
	bcw = bcleft.apply_left(w)
	Fx[ii,:bcw.shape[-1]] = bcw

#Right boundary conditionsbcright = BC_Mixed(X=x, bandwidth=bcbw, condition=m2, xbc=1)
for ii in range(Nx-bcend,Nx):
	xiis = ii-offset+bcorder; xiie = Nx+2*bcorder
	
	#Apply boundary conditions to stencil based on x with boundary
	#nodes included, but centered at ii+bcoffset
	w = fd_calc_weights(x[xiis:xiie],xint[ii],Nd)[Nd]
	bcw = bcright.apply_right(w)
	Fx[ii,-bcw.shape[-1]:] = bcw
tick.stop("sparse, BC")

Fx = Fx.tocsc()

tick.start("sparse, matvec")
#Apply Fx to the the correct dimension of v
y = Fx.matvec(v)
yr = Fx.rmatvec(v)
tick.stop("sparse, matvec")

Fas = LineSparse((Nx,Nx))

tick.start("alg, main")
#Difference matrix non-border elements
equalspaced=1
if equalspaced:
	w = fd_calc_weights(xint[:bw],xint[offset],Nd)[Nd]
	Fas.lines.default = Line(None,None,w,offset)
else:
	for ii in range(bcstart,Nx-bcend):
		iis = max(0,ii-offset); iie = min(ii+bw-offset,Nx)
		w = fd_calc_weights(xint[iis:iie],xint[ii],Nd)[Nd]
		Fas.lines.append(Line(ii,ii, w, offset))
tick.stop("alg, main")

tick.start("alg, BC")
#Left boundary conditions
for ii in range(0,bcstart):
	xiis = 0; xiie = ii+bw-offset+bcorder
	
	#Apply boundary conditions to stencil based on x with boundary
	#nodes included, but centered at ii+bcoffset
	w = fd_calc_weights(x[xiis:xiie],xint[ii],Nd)[Nd]
	bcw = bcleft.apply_left(w)
	Fas.lines.append(Line(ii,ii,bcw,ii))

#Right boundary conditions
for ii in range(Nx-bcend,Nx):
	xiis = ii-offset+bcorder; xiie = Nx+2*bcorder
	
	#Apply boundary conditions to stencil based on x with boundary
	#nodes included, but centered at ii+bcoffset
	w = fd_calc_weights(x[xiis:xiie],xint[ii],Nd)[Nd]
	bcw = bcright.apply_right(w)
	Fas.lines.append(Line(ii,ii,bcw,offset))
tick.stop("alg, BC")

tick.start("alg, matvec")
y2=Fas.matvec(v)
yr2=Fas.rmatvec(v)
tick.stop("alg, matvec")

#Test difference matrix
#w,v = linalg.eig(Fx.toarray())
w,v = sparse.linalg.eigen(Fx, k=5, which='SM')
print sqrt(-w)*X/pi

sigma = sqrt(-w)
print "DtN Error:",vstack(exp(2j*sigma*X)*(1j*sigma-m2)/(1j*sigma+m2)-(1j*sigma-m1)/(1j*sigma+m1))

import pylab

print "Check exact solution"
for kk in range(5):
	sigma = sqrt(-w[kk])
	xana = (sigma-1j*m1)*exp(1j*sigma*xint) + (sigma+1j*m1)*exp(-1j*sigma*xint)
	
	ev=median(v[:,kk]/xana)
	print abs(v[:,kk]/ev-xana).max() / abs(xana).max()

	pylab.plot(xint,real(v[:,kk]/ev),'r', xint,real(xana),'b--')

print tick.report()

