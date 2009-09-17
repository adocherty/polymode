#Unit test for hankel

from __future__ import division

import sys,os
sys.path.append(os.pardir)

import pylab as pl
from numpy import *
from numpy.testing import *

import timer
from scipy import random
from scipy.special import jv, hankel1, gamma

from PolyMode.mathlink.amos import hankel1_ratio_2d,  hankel1_ratio_1d, hankel1_ratio

def hankel1p(m, x):
	return 0.5*(hankel1(m-1,x)-hankel1(m+1,x))


def hankel_ratio_amos(m, x):
	"Calculates the ratio of Hankel functions: hankel_ratio = hankel1'(m,x)/hankel1(m,x)"
	#br = hankel1(m-1,x)/hankel1(m,x) - m/x
	if isscalar(m):
		br = hankel1_ratio(m-1, x)
	elif m.ndim==1:
		br = zeros(m.shape, dtype=complex128)
		hankel1_ratio_1d(m.astype(float64), x, br)
	elif m.ndim==2:
		br = zeros(m.shape, dtype=complex128)
		hankel1_ratio_2d(m.astype(float64), x, br)
	return br

def hankel_ratio_scipy(m, x):
	br = x*hankel1(m-1,x)/hankel1(m,x) - m
	return br

def hankel_ratio_scipy_p(m, x):
	br = x*hankel1p(m,x)/hankel1(m,x)
	return br

t = timer.timer()

Nt = 100
mr = random.random_integers(-500,500,(100,2))
xr = random.randn(Nt)*300 + random.randn(Nt)*30j

#"Exact" defined by scipy
exact_y = [hankel_ratio_amos(mr, x) for x in xr]

t.start("asym")
y = zeros(mr.shape, complex_)
for ii in range(Nt):
	z=xr[ii]
	m=mr
	d = (e*z/2/(m+1))**2*((m-1)/(m+1))**(m-0.5)
	b1 = m*(1-2*d)
t.stop("asym")

nans_found = error = 0
t.start("scipy")
for ii in range(Nt):
	y = hankel_ratio_scipy(mr, xr[ii])
	error = max(abs(y-exact_y[ii]).max(), error)
	nans_found = nans_found + any(isnan(y))
t.stop("scipy")

print "Scipy -- error:", error, "nans found:", nans_found

nans_found = error = 0
t.start("amos")
for ii in range(Nt):
	y = hankel_ratio_amos(mr, xr[ii])
	error = max(abs(y-exact_y[ii]).max(), error)
	nans_found = nans_found + any(isnan(y))
t.stop("amos")

print "Amos -- error:", error, "nans found:", nans_found

def h1_ratio(v,z):
	h1 = z*hankel1(v-1,z)/hankel1(v,z) - v
	return h1

#Asymptotic form for large v
def h1_ratio_asym(v,z,switch):
	m=abs(v)
	if switch==1:
		h1 = z**2/(2*m)-m
	elif switch==2:
		h1f=lambda m,z: 1/sqrt(2*pi*m)*(e*z/2/m)**m - 1j*sqrt(2)/sqrt(pi*m)*(2*m/e/z)**m 
		h1 = m*(h1f(m-1,z)-h1f(m+1,z))/(h1f(m-1,z)+h1f(m+1,z))
	elif switch==3:
		h1 = -m*(1-2*(z/2/(m+1))**2*(1+2/m+4/m**2/3)*sqrt((m-1)/(m+1)))
	else:
		d = sqrt((m+1)/(m-1))*(z/2)**2/(m+1)**2*(1+2/m+4/m**2/3)
		h1 = m*(d-1)/(d+1)
	return h1

def b1_ratio(v,z):
	h1 = z*jv(v-1,z)/jv(v,z) - v
	return h1

#Asymptotic form for large v
def b1_ratio_asym(v,z,switch):
	m=abs(v)
	if switch==1:
		b1 = m*((m+1)*(m+2)-(z/2)**2)/((m+1)*(m+2)+(z/2)**2)
	elif switch==2:
		b1 = m*(1-2*(z/2)**2/(m+1)/(m+2))
	elif switch==2:
		b1 = m*(sqrt(2*pi*(m+1))*(e*z/2/(m-1))**(m-1) - sqrt(2*pi*(m-1))*(e*z/2/(m+1))**(m+1)) /\
			(sqrt(2*pi*(m+1))*(e*z/2/(m-1))**(m-1) + sqrt(2*pi*(m-1))*(e*z/2/(m+1))**(m+1))
	elif switch==3:
		d = (e*z/2/(m+1))**2*((m-1)/(m+1))**(m-0.5)
		b1 = m*(1-2*d)
	elif switch==4:
		b1 = m*(1-2*(z/2/(m+1))**2*(1+2/m+4/m**2/3)*sqrt((m-1)/(m+1)))
	return b1

v = arange(-100,100)
z = 10+1j

pl.clf()

x = v
for ii in range(1,5):
	pl.plot(x,abs(h1_ratio(v,z) - h1_ratio_asym(v,z,ii))/abs(h1_ratio(v,z)),'--', label="H%d"%ii)
	pl.plot(x,abs(b1_ratio(v,z) - b1_ratio_asym(v,z,ii))/abs(b1_ratio(v,z)),'-', label="J%d"%ii)

pl.semilogy()
pl.legend()
pl.show()

print t.report()

