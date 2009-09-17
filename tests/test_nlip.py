from __future__ import division

import sys,os
if os.path.pardir not in sys.path:
	sys.path = [os.path.pardir] + sys.path

import logging
logging.basicConfig(level=logging.INFO, format='%(levelname).1s: %(message)s')

from pylab import *
from numpy import *
from scipy import linalg

from PolyMode import *
from PolyMode import KrylovSolver
reload(KrylovSolver)

#Solver.__polymode_debug__ = True

air = Material.Air()
mat = Material.Fixed(1.45)

wg = Waveguide.Waveguide(rmin=None, rmax=None, material=mat, symmetry=6)
c = Waveguide.Circle(air, radius=2.0)
wg.shapes = Waveguide.create_hexagonal_tiling(c, layers=1, D=6.75, symmetry=6)

Nx = 200,20						#number of radial & azimuthal points
wl = 1.0							#Wavelength
k0 = 2*pi/wl
m0 = 1
Nmodes = 5

all_convergence = []
all_modes = []

ksolver = KrylovSolver.NLBiArnoldiSolver(wg, Nx)
ksolver.initialize(wl, m0)
ksolver.create(generate=False)

if 'md1' not in locals():
	md1,md2 = ksolver(wl, m0, nefflist=[1.448,1.444])

eq = ksolver.equation
jac = ksolver.jacobian

mu1=79.9
mu2=79.5+0.1j
mu = 68+0.1j
v = random.random(prod(ksolver.base_shape+(2,))).astype(complex_)
u = random.random(prod(ksolver.base_shape+(2,))).astype(complex_)
if 0:
	eq.set_lambda(mu1, onlyev=0, onlyconst=0)
	eq.diff.generate()
	av=eq.matvec(v)

	eq.set_lambda(mu1, onlyev=0, onlyconst=1)
	eq.diff.generate()
	cv=eq.matvec(v)

	eq.set_lambda(mu1, onlyev=1, onlyconst=0)
	eq.diff.generate(leftbc=1,rightbc=1,nodefault=1)
	ev = eq.matvec(v)

	print "Total error:", abs(av-cv-ev).max()

def nlip(u,v,mu1,mu2):
	if abs(mu1-mu2)<1e-12:
		jac.set_lambda(mu1)
		Jv = jac.matvec(v)
		return dot(conj(u),v) - dot(conj(u),Jv)
	else:
		eq.set_lambda(mu1)
		Av = eq.matvec(v)
	
		eq.set_lambda(mu2)
		Av -= eq.matvec(v)
		return dot(conj(u),v) - dot(conj(u),Av)/(mu1-mu2)

def nlip2(u,v,mu1,mu2):
	if abs(mu1-mu2)<1e-12:
		jac.set_lambda(mu1)
		Jv = jac.matvec(v)
		return dot(conj(u),v) - dot(conj(u),Jv)
	else:
		eq.set_lambda(mu1, onlyev=1)
		Bv = eq.matvec(v)
	
		eq.set_lambda(mu2, onlyev=1)
		Bv -= eq.matvec(v)
		return dot(conj(u),v) - dot(conj(u),Bv)/(mu1-mu2)
	
eq.set_lambda(md1.evalue)
print "Residue:", abs(eq.matvec(md1.right)-md1.evalue*md1.right).max()
print "Left Residue:", abs(eq.rmatvec(md1.left)-conj(md1.evalue)*md1.left).max()

for md in [md1,md2]:
	print "Mode R NLIP:", nlip(md.left,v,md1.evalue,mu)
	print "Mode L NLIP:", nlip(u,md.right,mu,md1.evalue)

#Orthogonalize
ksolver.lorthogonalize = c_[md1.left, md2.left].T
ksolver.rorthogonalize = c_[md1.right, md2.right].T
ksolver.Northogonalize = 2

ksolver.approximate_orthogonalize(u,v,mu)
#ksolver.construct_alpha_vectors(u, v, al, ar)

for md in [md1,md2]:
	print "Mode R NLIP:", nlip(md.left,v,md.evalue,mu)
	print "Mode L NLIP:", nlip(u,md.right,mu,md.evalue)

