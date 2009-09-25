import pickle, os, sys
sys.path.append(os.path.pardir)

from numpy import *
from scipy.linalg import hessenberg
import householder
reload(householder)

from householder import *

def iprod(x,y):
	return dot(conj(x),y)
def inorm(x):
	return sqrt(iprod(x,x))

N=5
if 'V' not in locals():
	V = random.random( (N,N) )# + random.random( (N,N) )*1j
	#V[2] = V[1]+1e-12*V[2]

qcheck, rcheck = linalg.qr(V.T)
print "Residual for true:", absolute(dot(qcheck,rcheck)-V.T).max()

R = zeros((N,N), dtype=V.dtype)
hhqr = HouseholderQR(overwrite = True)
q1 = zeros(N, dtype=V.dtype)

Vx = V.copy()
hhqr(Vx, R)

Q = eye(N, dtype=V.dtype)
hhqr.q_action( Q )

print "HH QR Residual:", absolute(dot(Q,R)-V.T).max()
print "Orthogonality:", absolute(iprod(Q.T,Q)-eye(N)).max()

#absolute(R-V.T).max()

QV=V.copy()
hhqr.q_action(QV)
print "Difference in Qx", absolute( dot(Q,V) - QV ).max()

QV=V.copy()
hhqr.q_inverse(QV)
print "Difference in Q*x", absolute( dot(conj(Q.T),V) - QV ).max()

QV=V.copy()
hhqr.q_raction(QV)
print "Difference in xQ", absolute( dot(V,Q) - QV ).max()

QV=V.copy()
hhqr.q_rinverse(QV)
print "Difference in xQ*", absolute( dot(V,conj(Q.T)) - QV ).max()

#Check Householder with zero row
Vx = V.copy()
Vx[0,0:]=0

hhqr = HouseholderQR(overwrite = False)
hhqr(Vx,R)
Q = eye(N, dtype=V.dtype); hhqr.q_action( Q )
print "Householder on zero matrix:", absolute(dot(Q,R)-Vx.T).max()

#Check Householder with zero row
Vx = V.copy()
Vx[0:,0]=0

hhqr = HouseholderQR(overwrite = False)
hhqr(Vx,R)
Q = eye(N, dtype=V.dtype); hhqr.q_action( Q )
print "Householder on zero pivot:", absolute(dot(Q,R)-Vx.T).max()

#Check Householder with zero row
Vx = V.copy()
Vx[0,0:]=0

hhqr = HouseholderQR(overwrite = False)
hhqr(Vx,R)
Q = eye(N, dtype=V.dtype); hhqr.q_action( Q )
print "Householder on zero row:", absolute(dot(Q,R)-Vx.T).max()


#Try Hessenburg reductions
A = random.random((N,N))
A[:,0]=0
A[0,:]=0

hhh=HouseholderHessenberg(overwrite=True)
Ax = A.copy()
H = hhh(Ax)
print H

Q = eye(N, dtype=complex_)
hhh.q_action(Q)

A2 = A.copy()
hhh.q_raction( A2 )
hhh.q_inverse( A2 )
print "Hessenburg Action", absolute(A2 - H).max()

hrow = RowHessenberg(overwrite=True)
Ax = A.copy()
Hx = hrow(Ax)
Qx = eye(N, dtype=complex_)
hrow.q_action(Qx)
print "Hessenburg Row Action", absolute(dot(conj(Qx.T), dot(A,Qx)) - Hx).max()

#Reverse hessenberg transforms
hrow.q_action(Hx.T)
hrow.q_rinverse(Hx.T)
print "Hessenburg Row Reverse Action", absolute(A - Hx).max()

w = zeros(N, dtype=complex_); w[-1]=1
hrow.q_raction(w)
hrow.q_rinverse(w)
hrow.q_action(w[:,newaxis])
hrow.q_inverse(w[:,newaxis])
print "Row final transform",w

if 0:
	ghr=HessenbergReduce(overwrite=True)
	ghr(H)
	ghr.q_action( A2 )
	print "Hessenburg Reduction Check" , absolute(A2-H).max()


	#Test deflation of an eigenvector of X so W* X W e1 = l e1
	# ie the eigenvector is transformed to be e1

	X = random.random((N,N))+random.random((N,N))*1j
	evals, evecs = linalg.eig(X)

	print "Deflating ",evals[0:3]
	hh = HouseholderQR(overwrite=False)
	hh(evecs[:,:3].T)

	hh.q_raction(X)
	hh.q_inverse(X)
	print "Errors in deflated evalues", absolute(sort(diag(X)[:3])-sort(evals[:3])).max()


