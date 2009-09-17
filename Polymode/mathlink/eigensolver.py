# coding=utf-8
'''
┈┈┈┈┈┈┈┈┈┈ Suite for solution of the linear eigenproblem  ┈┈┈┈┈┈┈┈┈┈
A python/numpy implementation of the Implicitly Restarted Arnoldi Method (IRAM)

Still to do: locking/purging of converged eigenvectors

┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈
'''

from __future__ import division

from numpy import *
import os, sys, logging

from misc import machine_precision
from householder import *

def iprod(x,y):
	return dot(conj(x),y)
def inorm(x):
	return sqrt(iprod(x,x))

#
# Some paramters for the eigensolver
#
eigensolver_debug = False
eigensolver_deflation = True

#
# Arnoldi solve must be given a class with matvec and rmatvec members
# These are some standard such matrix classes
#
class ArnoldiDense(object):
	def setup(self, Ain):
		self.A = Ain
		self.shape = Ain.shape
		self.dtype = Ain.dtype

	def residue(self, v, u):
		return absolute(dot(self.A, v) - u*v).max()

	def matvec(self, x):
		y = dot(self.A, x)
		return y

	def rmatvec(self, x):
		y = dot(x.T, self.A)
		return y

class ArnoldiShiftInvert(object):
	def eigenvalue_transform(self, hvals):
		"Map the shift-invert eigenvalues to the eigenvalues of A"
		return 1./hvals + self.shift

	def residue(self, v, u):
		return absolute(dot(self.A, v) - self.eigenvalue_transform(u)*v).max()

	def matvec(self, x):
		return x

	def rmatvec(self, x):
		return x

class ArnoldiShiftInvertDense(ArnoldiShiftInvert):
	def __init__(self, overwrite=False):
		#Import linalg routines from scipy
		from scipy import linalg
		self.linalg = linalg
		self.overwrite = overwrite

	def setup(self, Ain, shift=0.0):
		"Setup LU factorization A-sI = LU with a shift s"
		self.shift = shift
		self.A = Ain
		self.shape = Ain.shape[0],Ain.shape[0]
		self.dtype = Ain.dtype
		self.lu = self.linalg.lu_factor(Ain - shift*eye(Ain.shape[0]))

	def residue(self, v, u):
		return absolute(dot(self.A, v) - self.eigenvalue_transform(u)*v).max()

	def matvec(self, x):
		y = self.linalg.lu_solve(self.lu, x)
		return y

	def rmatvec(self, x):
		y = self.linalg.lu_solve(self.lu, x, trans=1)
		return y

class ArnoldiShiftInvertGmres(ArnoldiShiftInvert):
	def __init__(self):
		from scipy import linalg
		self.linalg = linalg

	#For shift & invert A must be CSC/CSR
	def setup(self, matrix, shift=0.0):
		self.matrix = matrix
		self.shift = shift
		self.shape = matrix.shape
		self.dtype = matrix.dtype

	def matvec(self, x):
		xn = linalg.gmres(self.matrix, x, tol=1e-10)[0]
		return xn

	def rmatvec(self, x):
		pass

#
# To choose which eigenvalues are wanted and which are unwanted
# to be applied to the implict restart routine, use one of these
# classes
#
class EigenvalueSelector(object):
	pass
	
class ClosestAbsolute(EigenvalueSelector):
	'''
	Selector to choose eigenvalues closest to the shift
	(i.e. largest magnitude)
	'''
	def __call__(self,shifts,P):
		isort = absolute(shifts).argsort()
		iunwanted, iwanted = split(isort,[-P])
		return iunwanted, iwanted

class ClosestReal(EigenvalueSelector):
	'''
	Selector to choose eigenvalues with real part closest from the shift
	(i.e. largest real part)
	'''
	def __call__(self,shifts,P):
		isort = absolute(real(shifts)).argsort()
		iunwanted, iwanted = split(isort,[-P])
		return iunwanted, iwanted

class FarthestAbsolute(EigenvalueSelector):
	'''
	Selector to choose eigenvalues farthest from the shift
	(i.e. smallest magnitude)
	'''
	def __call__(self,shifts,P):
		isort = absolute(shifts).argsort()
		iwanted, iunwanted = split(isort,[P])
		return iunwanted, iwanted

class FarthestReal(EigenvalueSelector):
	'''
	Selector to choose eigenvalues with real part farthest
	from the shift (i.e. smallest real part)
	'''
	def __call__(self,shifts,P):
		isort = absolute(real(shifts)).argsort()
		iwanted, iunwanted = split(isort,[P])
		return iunwanted, iwanted

#For convenience define other names
LargestMagnitude = ClosestAbsolute
LargestReal = ClosestReal
SmallestMagnitude = FarthestAbsolute
SmallestReal = FarthestReal

# Orthonormalizing routines that construct column jj of H
# by orthogonalizing column jj+1 against columns 0..jj of V
def mgs(V, H, jj):
	"""
	 Perform single orthogonalization step using Modified Gram-Schmidt with
	 
	 Inputs:
	  V: Vectors
	  H: Dense storage of Hessian projection matrix
	  jj: orthogonalize jth vector against vectors 0..j-1
	"""
	w = V[jj+1]
	for ii in range(0,jj+1):
		h = iprod(V[ii],w)
		w = w - h*V[ii]
		H[ii,jj] = h

	#Iterate reorthogonalization until orthogonality is acheived
	hmax=abs(h); ni=0
	while hmax>1e-10 and ni<5:
		ni+=1
		hmax=0
		for ii in range(0,jj+1):
			h = iprod(V[ii],w)
			w = w - h*V[ii]
			H[ii,jj] += h
			hmax=max(abs(h),hmax)

	hjp = inorm(w)
	V[jj+1] = w/hjp

	return hjp

def djks(V, H, jj):
	"""
	 Perform single orthogonalization step using Gram-Schmidt with
	 Daniel, Gragg, Kaufman, and Stewart refinement
	 This implementation allows for BLAS level 2 speed-up
	 
	 Inputs:
	  V: Vectors
	  H: Dense storage of Hessian projection matrix
	  jj: orthogonalize jth vector against vectors 0..j-1
	"""

	w = V[jj+1]
	h = iprod(V[:jj+1], w)
	w = w - dot(h, V[:jj+1])

	#Reorthogonalize if needed
	eta = 1/sqrt(2)
	iterations = 0
	s = iprod(V[:jj+1], w)
	while (iterations<5) and absolute(s).max()>1e-9:
		iterations+=1
		w = w - dot(s,V[:jj+1])
		h = h + s
		s = iprod(V[:jj+1], w)

	H[:jj+1,jj]=h
	hjp = inorm(w)
	V[jj+1] = w/hjp
	return hjp	

class ArnoldiIR:
	'''
	Implicitly restarted Arnoldi algorithm
	'''
	def __init__(self, matrix, N, Nv, selector=ClosestAbsolute(), dtype=complex):
		self.matrix = matrix
		self.N = N
		self.selector = selector
		self.dtype = dtype
		
		##Set up V & H space
		self.Nv = Nv
		self.Vm = empty( (Nv+1,N), dtype=self.dtype )
		self.Hm = zeros( (Nv+1,Nv), dtype=self.dtype )
		
		self.tolerance = 1e-10
		self.maxit = 10
		self.transpose = False
		
		self.kmin = 0

	def rayleigh_ritz_factor(self,jstart,jend):
		"""
		Arnoldi factorization of A
		  A Vm = Vm Hm + fm em^T
		Start at column jstart and go to column jend
		"""
		hjp = 1
		#Normalize start vector
		vnrm = inorm(self.Vm[jstart])
		self.Vm[jstart] /= vnrm

		for jj in range(int(jstart),int(jend)):
			if self.transpose==1:
				self.Vm[jj+1] = self.matrix.rmatvec(self.Vm[jj])
			elif self.transpose==2:
				self.Vm[jj+1] = conj(self.matrix.rmatvec(conj(self.Vm[jj])))
			else:
				self.Vm[jj+1] = self.matrix.matvec(self.Vm[jj])

			hjp = mgs(self.Vm, self.Hm, jj)
			self.Hm[jj+1,jj]=hjp

		self.hjp_last = hjp

	def Heig(self,n=None):
		if n:
			self.Hval,self.Hvec = linalg.eig(self.Hm[:n,:n])
		else:
			self.Hval,self.Hvec = linalg.eig(self.Hm[:-1, :])
		return self.Hval

	def ritz_values(self,iselect):
		return self.Hval[iselect]

	def ritz_vectors(self, iselect):
		ar_vec = dot(self.Hvec[:,iselect].T, self.Vm[:-1])
		return ar_vec

	def residue(self, iselect):
		return absolute(self.Hvec[-1,iselect]*self.hjp_last)

	#Implicit restart
	def restart_old(self, ushifts):
		kmin = self.kmin

		#Length of restarted Arnoldi problem:
		P = self.Nv - len(ushifts)

		Hm = self.Hm[kmin:-1,kmin:]
		Im = eye(self.Nv-kmin, dtype=self.dtype)
		q = zeros(self.Nv-kmin, dtype=self.dtype); q[-1]=1

		#Implicit QR shifts on Hm for unwanted shifts
		for mu in ushifts:
			Qj,Rj = linalg.qr(Hm - mu*Im)
			
			Hm[:] = dot(conj(Qj).T, dot(Hm,Qj))
			self.Hm[:kmin, kmin:] = dot(self.Hm[:kmin, kmin:],Qj)

			self.Vm[kmin:-1] = dot(Qj.T, self.Vm[kmin:-1])
			q = dot(Qj.T, q)

		#Update residual vector:
		fm = self.Vm[P]*self.Hm[P,P-1] + self.Vm[-1]*self.Hm[-1,-1]*q[P-kmin-1]
		self.Hm[P,P-1] = inorm(fm)
		self.Vm[P] = fm/self.Hm[P,P-1] 

		#Fix up zero Hm entries
		self.Hm[P+1:]= 0
		return P


	#Implicit restart
	def restart(self, ushifts):
		kmin = self.kmin

		#Length of restarted Arnoldi problem:
		P = self.Nv - len(ushifts)
		Im = eye(self.Nv-kmin, dtype=self.dtype)
		q = zeros(self.Nv, dtype=self.dtype); q[-1]=1

		#Implicit QR shifts on Hm for unwanted shifts
		hh = HouseholderQR(overwrite=False)
		for mu in ushifts:
			hh.clear_transforms()

			#Compute orthogonal factorization of Xi = WR
			hh(self.Hm[kmin:-1, kmin:].T - mu*Im)

			#Update	H+ = Qmu* Hm Qmu
			#Note Hm-j+ = Qmu* Hm-j Qmu and Gm-j+ = Gm-j Qmu
			#The action does this implicitly
			hh.q_inverse(self.Hm[:-1, :])
			hh.q_raction(self.Hm[:-1])

			#Vm-j+ = Vm-j Q
			hh.q_raction(self.Vm[:-1].T)
			hh.q_raction(q)

		#Update residual vector fk = vk Hm(k,k-1) + fm em^T Q
		fm = self.Vm[P]*self.Hm[P,P-1] + self.Vm[-1]*self.hjp_last*q[P-1]
		self.Hm[P,P-1] = inorm(fm)
		self.Vm[P] = fm/self.Hm[P,P-1] 

		#Fix up zero Hm entries
		self.Hm[P+1:]= 0
		return P
		
	#Detect any vectors that have low residue and lock them
	def deflate(self):
		res = absolute(self.Hvec[-1, self.kmin:]*self.hjp_last)
		
		#Select vectors to deflate
		ilock = nonzero(res<self.tolerance)[0]
		
		#There is a strange error that causes the eigensolver to deflate the last
		#entry after previous deflations, why does this have such small residue?
		if len(ilock)>0 and ilock[0]<len(ilock):
			if eigensolver_debug: print "Deflating %s" % (ilock,)
			Xi = self.Hvec[self.kmin:, self.kmin:][:, ilock]
			self.lock(Xi)
			self.kmin += len(ilock)
	
	#lock the eigenvectors 
	def lock(self, Xi):
		i = shape(Xi)[1]
		kmin = self.kmin

		#Compute orthogonal factorization of Xi = WR
		qr = HouseholderQR(overwrite=True)
		qr(Xi.T)

		#Update	H+ = W* Hm W
		qr.q_inverse(self.Hm[:-1])
		qr.q_raction(self.Hm[:-1])

		#Update	Vm+ = Vm W
		qr.q_raction(self.Vm[kmin:-1].T)

		#Compute P to restore Hm to hessenburg form
		hh = RowHessenberg(overwrite=True)
		hh(self.Hm[kmin+i:-1,kmin+i:])
		
		#Update masked part of Hm++ = Hm+ Qh
		hh.q_inverse(self.Hm[:kmin+i, kmin+i:].T)
		
		#Update	Vm++ = Vm+ Qh
		hh.q_inverse(self.Vm[kmin+i:-1])
		
	def check(self, which, unwhich):
		Hm=self.Hm
		Vm=self.Vm
		
		print
		
		#Test esolution
		avec = self.ritz_vectors(which)
		val = self.ritz_values(which)
		fres_max = 0; fres_min=inf
		for ii in range(len(which)):
			fres = absolute(self.matrix.matvec(avec[ii]) -  val[ii]*avec[ii]).max()
			fres_min = min(fres_min,fres); fres_max = max(fres_max,fres)
		print "Wanted Residue: %.4g -> %.4g" % (fres_min, fres_max)

		uvec = self.ritz_vectors(unwhich)
		uval = self.Hval[unwhich]
		ures_max = 0; ures_min=inf
		for ii in range(len(unwhich)):
			ures = absolute(self.matrix.matvec(uvec[ii]) -  uval[ii]*uvec[ii]).max()
			ures_min = min(ures_min,ures); ures_max = max(ures_max,ures)
		print "Unwanted Residue: %.4g -> %.4g" % (ures_min, ures_max)

		#Check arnoldi expansion with residue
		X = self.Vm[:-1].T.copy()
		for i in range(X.shape[1]):
			X[:,i] = self.matrix.matvec(X[:,i])
		X -= dot(self.Vm[:-1].T, self.Hm[:-1])
		X[:,-1] -= self.Hm[-1,-1]*self.Vm[-1]
		print "Error in Arnoldi expansion:", absolute(X).max()

		print "V orthogonality: %.4g" % absolute(iprod(Vm,Vm.transpose())-eye(self.Nv+1)).max()

		#print "Res:",absolute(self.Hvec[-1]*self.hjp_last)
		
		return fres

	def __call__(self, number, v0=None):
		"""
		Run the IRAM to calculated number largest eigenvalues using v0 as a starting vector
		"""
		#Number of eigenvalues to find
		P = number
		M = self.Nv

		self.tracked_solutions = []
		
		#Starting vector
		if v0==None:
			self.Vm[0] = (random.rand(self.N)+random.rand(self.N)*1j).astype(self.dtype)
		else:
			self.Vm[0] = v0

		self.rayleigh_ritz_factor(0, M)					#Generate initial Arnoldi basis
		
		ar_conv = inf; niter=0
		while self.kmin<P and ar_conv>self.tolerance and niter<self.maxit:
			#Calculate eigeninformation
			Hval = self.Heig()
			iunwanted, iwanted = self.selector(Hval, P)

			if eigensolver_debug: self.check(iwanted, iunwanted)

			#Convergence estimates
			ar_residue = self.residue(iwanted)
			ar_conv = ar_residue.max()
	
			if len(iunwanted)>0 and ar_conv>self.tolerance:
				if eigensolver_deflation:
					self.deflate()					#Deflation of converged vectors
				P = self.restart(Hval[iunwanted])	#Implicit restart using exact shifts
				self.rayleigh_ritz_factor(P, M)		#Expand Arnoldi basis

			niter+=1
			if eigensolver_debug:
				print "[%d] Res: %.4g->%.4g, Deflated: %d" \
					% (niter, ar_residue.max(), ar_residue.min(), self.kmin)

			#Save all converging solutions for debugging
			self.tracked_solutions += [ (self.ritz_values(iwanted), \
				self.ritz_vectors(iwanted)) ]

		#Ritz vectors and values
		ar_eval = self.ritz_values(iwanted)
		ar_vec = self.ritz_vectors(iwanted)
		self.final_residue = ar_residue
		
		logging.debug( "IRAM: Converged to: %.4g in %d iterations" % (ar_conv, niter) )

		return ar_eval, ar_vec

class HouseholderArnoldiIR(ArnoldiIR):
	"""
	Implicitly Restarted Arnoldi Method
	implemented with Householder transformations
	"""
	def __init__(self, *args, **kwargs):
		ArnoldiIR.__init__(self, *args, **kwargs)
		self.irtransforms = transform_action()
		
	def rayleigh_ritz_factor(self,kstart,kend):
		#Careful of the order in whcih the Q's are multiplied!!
		orthogonalize = householder_qr(overwrite=True, offset=1)
		
		for jj in range(kstart,kend+1):
			#Apply prior IR transforms to vectors
			Vtmp = self.Vm[jj:jj+1]
			self.irtransforms.q_inverse( Vtmp.T )
			orthogonalize.q_inverse( Vtmp.T )
			
			#Construct column jj of Hm
			orthogonalize.qr(self.Vm, self.Hm, jj, jj+1)
			
			#Generate new orthogonal basis vector jj
			self.Vm[jj] = 0; self.Vm[jj,jj]=1
			orthogonalize.q_action( self.Vm[jj] )
			self.irtransforms.q_action( self.Vm[jj] )

			#Construct next vector by action
			if jj<self.Nv:
				self.Vm[jj+1] = self.matrix.matvec(self.Vm[jj])

		#Save to list of transforms
		self.irtransforms.add_transforms(orthogonalize.transforms)

		self.hjp_last = self.Hm[-1,-1]

	def generate_V(self, jj):
		#Generate new orthogonal basis vector jj
		self.Vm[jj] = 0; self.Vm[jj,jj]=1
		self.irtransforms.q_action( self.Vm[jj] )

	def restart(self, ushifts, P):
		kmin = self.kmin
		Nvm = self.Nv-kmin

		Hm = self.Hm[kmin:-1,kmin:]
		Im = eye(Nvm, dtype=self.dtype)

		irqr = householder_qr()
		Qir = transform_action()
		
		#Implicit QR shifts on Hm for unwanted shifts
		for mu in ushifts:
			irqr.qr(Hm - mu*Im)
			irqr.q_inverse( Hm )				#Update Hm with action of Q
			irqr.q_raction( Hm )

			#We'll add the IR transforms to this group
			Qir.add_transforms( irqr.transforms, offset=Nvm )
			irqr.transforms = []

		#Calculate element of Q needed for residual
		q = zeros(Nvm, dtype=self.dtype); q[P-1]=1
		Qir.q_action( q )

		#Update residual vector:
		fm = self.Vm[P]*self.Hm[P,P-1] + self.Vm[-1]*self.hjp_last*q[-1]
		self.Hm[P,P-1] = inorm(fm)
		self.Vm[P] = fm/self.Hm[P,P-1] 
	
		#Add QR decomposition to main transforms
		self.irtransforms.add_transforms( Qir.transforms )

		#Fix up zero Hm entries
		self.Hm[P+1:]= 0


def eigs(A, number=None, shift=None, v0=None, tol=1e-10, maxiter=100, \
		selector=ClosestAbsolute(), transpose=0):
	"""
	Calculate a few eigenvalues of A using the IRAM
	A: matrix
	shift: use shift and invert to caclulate eigenvalues
	number: how many eigenpairs to generate
	tol: tolerance for convergence
	transpose: solve right eigenproblem
	
	returns: (u, v)
	u: eigenvalues
	v: eigenvectors
	"""
	#if number is not specified make it 10
	if number is None:
		number = min(A.shape[0]/2, 10)
	
	#Guess how many Arnoldi basis vectors to use
	M = min(max(10, 3*number), A.shape[0])
	
	#Use A as matrix class if it has a matvec member
	if hasattr(A, 'matvec'):	
		matrix = A

	#Use dense Arnoldi class (assume dense matrix)
	elif shift is None:
		matrix = ArnoldiDense()
		matrix.setup(A)

	#Use shift invert if shift given (assume dense matrix)
	else:
		matrix = ArnoldiShiftInvertDense(overwrite=False)
		matrix.setup(A, shift=shift)

	#Instanciate a solver and run it
	#selector = ClosestAbsolute()
	asolve = ArnoldiIR(matrix, A.shape[0], M, selector=selector, dtype=A.dtype)
	asolve.transpose = transpose
	asolve.maxit = maxiter
	asolve.tol = tol
	evals, evecs = asolve(number, v0=v0)
	
	#Get the corrected eigenvalues
	if hasattr(matrix, 'eigenvalue_transform'):
		evals = matrix.eigenvalue_transform(evals)

	if asolve.final_residue.max()>tol:
		print "Warning: Eigensolve has not converged to required tolerance"
		print asolve.final_residue

	#Return the calculated eigenvalues and vectors
	return evals, evecs
	
	
