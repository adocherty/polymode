"""
Bugs:
* 
"""


import logging, pickle, os, sys
logging.basicConfig(level=logging.DEBUG, format='%(levelname).1s: %(message)s')
sys.path.append(os.path.pardir)

from numpy import *
import eigensolver, householder
reload(eigensolver)

set_printoptions(precision=4, linewidth=100, suppress=0)

#Test eigensolver dense matrix solve
Na = 1000
neigs = 20
datatype = complex128

if "A" not in locals():
	A = random.random((Na,Na)) + 0*1j
	v0 = random.random(Na) + 0*1j

shift = random.rand() + random.rand()*1j
evals, evecs = eigensolver.eigs(A, neigs, shift=shift, v0=v0, maxiter=100)

#Check residue
for l,v in zip(evals,evecs):
	print "V Norm: %.2g,  Final Residue: %.2g" %( sqrt(dot(conj(v),v)), absolute(dot(A,v) - l*v).max() )

