# _*_ coding=utf-8 _*_
'''
Suite of tools for differentiation matricies
for both dense and banded matricies

 Changes:

 To Do:
'''
from __future__ import division
from numpy import *
from .line_sparse import Line, LineMulti, LineSparse

def fd_calc_weights(x, x0, m, dtype=float_):
	'''
	Calculate the Finite Difference weights for a FD stencil on a base of points x
	at the locations x0. Calculates m levels of stencils (m=0 corresponding to the identity)
	
	Algorithm taken from:
	B. Fornberg, “Generation of Finite Difference Formulas on Arbitrarily Spaced Grids,”
	Mathematics of Computation 51, no. 184 (1988): 699-706.
	'''
	x=asanyarray(x)
	n = len(x)

	weights = zeros((m+1,n), dtype=dtype)
	weights[0,0] = 1.
	betaold = 1.
	for i in range(1,n):
		beta = prod(x[i]-x[0:i])
		for k in range(0,min(i,m)+1):
			weights[k,i] = betaold*(k*weights[k-1,i-1]-(x[i-1]-x0)*weights[k,i-1])/beta
		betaold=beta

		for j in range(0,i):
			for k in range(min(i,m),-1,-1):
				weights[k,j] = ((x[i]-x0)*weights[k,j]-k*weights[k-1,j])/(x[i]-x[j])

	# Clear very small entries:
	weights[absolute(weights)<1e-10] = 0
	return weights

class BoundaryCondition:
	bcorder=0
	row_extent=1
	def apply_left(self, w):
		return w
	def apply_right(self, w):
		return w

class FiniteDifference:
	'''
	if X=vector, taken as the xnodes
	If X=scalar, xnodes are constructed as N equally spaced nodes from 0 ... X-dx
	'''
	def __init__(self, dmax=2, bandwidth=3, X=1, N=1, 
						bcl=None, bcr=None, dtype=complex_):
		self.dmax = dmax
		self.bandwidth = bandwidth
		self.offset = bandwidth//2
		self.dtype = dtype
		
		#Boundary conditions - would be nice to set the X here not in the
		#instanciation of the BC
		self.bcleft = BoundaryCondition() if bcl is None else bcl
		self.bcright = BoundaryCondition() if bcr is None else bcr
		self.equispaced = 1
		
		#If X is a list then use it as the nodes, otherwise take x = 0..X
		if isinstance(X, (ndarray, list)):
			self.xnode = X
		else:
			dx = X/N
			self.xnode = arange(0,N)*dx
		
		self.Nx = len(self.xnode)-self.bcleft.bcorder-self.bcright.bcorder
		self.Nx_bc = len(self.xnode)

	def __str__(self):
		rs = "Finite Difference with bandwidth %d\n | Right BC: %s\n | Left BC: %s" \
			% (self.bandwidth, self.bcleft, self.bcright)
		return rs

	def bc_extents(self):
		bcstart = self.offset
		bcend = self.bandwidth-self.offset-1
		return (bcstart, bcend)

	def stencil(self, nodes, at):
		weights = fd_calc_weights( self.xnode[nodes], self.xnode[at], self.dmax, dtype=self.dtype)
		return weights

	def set_right_bc(self, bc, chain=1):
		self.bcright.set_condition(bc, chain)

	def set_left_bc(self, bc, chain=1):
		self.bcleft.set_condition(bc, chain)

	def diff_at(self, x, at, n=1):
		weights = fd_calc_weights( self.xnode, at, self.dmax, dtype=self.dtype)[n]
		return dot(x,weights)

	#Some quick functions
	def diff0(self, x):
		return self.diff(x,0)
	def diff1(self, x):
		return self.diff(x,1)
	def diff2(self, x):
		return self.diff(x,2)

	def rdiff0(self, x):
		return self.rdiff(x,0)
	def rdiff1(self, x):
		return self.rdiff(x,1)
	def rdiff2(self, x):
		return self.rdiff(x,2)

	def diff(self, f, n=1, axis=0):
		"Equispaced difference without boundary condtions"
		fd_ord=self.fd_ord
		weights = self.stencil(slice(None,fd_ord+1), fd_ord//2)[dn]
		df = correlate( f.ravel(), weights, mode=1 )
		return df

	def rdiff(self, f, n=1, axis=0):
		return self.diff(f,dn)

class SparseDifferenceMatrix(FiniteDifference):
	def __init__(self, dmax, **kwargs):
		FiniteDifference.__init__(self, **kwargs)
		self.dn = dmax
		from scipy import sparse
		self.matrix = sparse.lil_matrix((self.Nx,self.Nx), dtype=self.dtype)

		#Create the FD class storage
		self.generate()

	def generate(self, leftbc=0, rightbc=0):
		full = not leftbc and not rightbc

		Nx = self.Nx
		ox = self.offset
		bw = self.bandwidth
		bco = self.bcleft.bcorder
		
		#Extents over which the boundary condition changes the lines
		bcstart,bcend = self.bc_extents()

		if full:
			for ii in range(bcstart,Nx-bcend):
				iis = max(0,ii-ox); iie = min(ii+bw-ox,Nx)
				self.matrix[ii, iis:iie] = self.stencil(slice(iis+bco,iie+bco),ii+bco)[self.dn]

		if full or leftbc:
			#Left boundary conditions
			for ii in range(0,bcstart):
				xii = ii+bco
				xiis = 0; xiie = xii+bw-ox
	
				#Apply boundary conditions to stencil based on x with boundary
				#nodes included, but centered at ii+bcoffset
				w = self.stencil(slice(xiis,xiie),xii)[self.dn]
				bcw = self.bcleft.apply_left(w)
				self.matrix[ii,:bcw.shape[-1]] = bcw

		if full or rightbc:
			#Right boundary conditions
			for ii in range(Nx-bcend,Nx):
				xii = ii+bco
				xiis = xii-ox; xiie = self.Nx_bc
	
				#Apply boundary conditions to stencil based on x with boundary
				#nodes included, but centered at ii+bcoffset
				w = self.stencil(slice(xiis,xiie),xii)[self.dn]
				bcw = self.bcright.apply_right(w)
				self.matrix[ii,-bcw.shape[-1]:] = bcw
		
		return self.matrix

class DifferenceMatrix(FiniteDifference):
	"""
	Finite difference class for generating and applying finite difference stencils
	with boundary conditions.
	"""
	def __init__(self, dmax, **kwargs):
		FiniteDifference.__init__(self,**kwargs)
		
		self.dmatrix = {}
		for dn in range(self.dmax+1):
			self.dmatrix[dn] = LineSparse((self.Nx,self.Nx))

		#Create the FD class storage
		self.generate()

	def __getstate__(self):
		"Pickle all needed data, ignore cached data"
		state = self.__dict__.copy()
		state['dmatrix'] = {}
		return state
	
	def __setstate__(self,state):
		"Restore pickled data"
		self.__dict__.update(state)
		
		self.dmatrix = {}
		for dn in range(self.dmax+1):
			self.dmatrix[dn] = LineSparse((self.Nx,self.Nx))

	def generate(self, leftbc=0, rightbc=0, nodefault=0):
		"Generate the internal stencils for the finite difference matrix"
		full = not leftbc and not rightbc
		
		Nx = self.Nx
		ox = self.offset
		bw = self.bandwidth
		bco = self.bcleft.bcorder
		
		#The difference stencils to create
		dns = range(self.dmax+1)
		
		#Extents over which the boundary condition changes the lines
		bcstart,bcend = self.bc_extents()

		if full and self.equispaced:
			w = self.stencil(slice(0,bw),ox)
			for dn in dns:
				self.dmatrix[dn].lines.default = Line(w[dn],ox)

		elif full:
			for ii in range(bcstart,Nx-bcend):
				iis = max(0,ii-ox); iie = min(ii+bw-ox,Nx)
				w = self.stencil(slice(iis+bco,iie+bco),ii+bco)

				for dn in dns:
					self.dmatrix[dn].lines[ii] = Line(w[dn], ox)
		elif nodefault:
			for dn in dns:
				self.dmatrix[dn].lines.default = None


		if full or leftbc:
			#Left boundary conditions
			for ii in range(0,bcstart):
				xii = ii+bco
				xiis = 0; xiie = xii+bw-ox
	
				#Apply boundary conditions to stencil based on x with boundary
				#nodes included, but centered at ii+bcoffset
				w = self.stencil(slice(xiis,xiie),xii)
				for dn in dns:
					bcw = self.bcleft.apply_left(w[dn])
					if bcw.ndim==1:
						self.dmatrix[dn].lines[ii] = Line(bcw, ii)
					else:
						self.dmatrix[dn].lines[ii] = LineMulti(bcw, ii)

		if full or rightbc:
			#Right boundary conditions
			for ii in range(Nx-bcend,Nx):
				xii = ii+bco
				xiis = xii-ox; xiie = self.Nx_bc
	
				#Apply boundary conditions to stencil based on x with boundary
				#nodes included, but centered at ii+bcoffset
				w = self.stencil(slice(xiis,xiie),xii)
				for dn in dns:
					bcw = self.bcright.apply_right(w[dn])
					if bcw.ndim==1:
						self.dmatrix[dn].lines[ii] = Line(bcw, offsetr=Nx-ii)
					else:
						self.dmatrix[dn].lines[ii] = LineMulti(bcw, offsetr=Nx-ii)

	#Application of the differences
	def diff(self, x, n=1, axis=0):
		"Apply the difference matrix D_n to vector x"
		x = x.astype(complex_)
		return self.dmatrix[n].matvec(x.swapaxes(0,axis)).swapaxes(axis,0)

	def rdiff(self, x, n=1, axis=0):
		"Apply the adjoint difference matrix (D_n)^A to vector x"
		if n==0: return x
		return self.dmatrix[n].rmatvec(x.swapaxes(0,axis)).swapaxes(axis,0)

	#Used to construct a difference matrix using this class
	#Think of a nicer way...
	def construct(self, ii):
		class FakeFD(object):
			def __init__(self, D,ii):
				self.ii = ii
				self.dmatrix=D
			def diff0(self, x):
				return self.dmatrix[0].mirror(x,self.ii)
			def diff1(self, x):
				return self.dmatrix[1].mirror(x,self.ii)
			def diff2(self, x):
				return self.dmatrix[2].mirror(x,self.ii)
		return FakeFD(self.dmatrix,ii)

class JacobianMatrix(DifferenceMatrix):

	def generate(self, leftbc=0, rightbc=0):
		"Generate the internal stencils for the finite difference matrix"
		Nx = self.Nx
		ox = self.offset
		bw = self.bandwidth
		bco = self.bcleft.bcorder
		
		#The difference stencils to create
		dns = range(self.dmax+1)
		
		#Extents over which the boundary condition changes the lines
		bcstart,bcend = self.bc_extents()
		
		if leftbc:
			#Left boundary conditions
			for ii in range(0,bcstart):
				xii = ii+bco
				xiis = 0; xiie = xii+bw-ox
	
				#Apply boundary conditions to stencil based on x with boundary
				#nodes included, but centered at ii+bcoffset
				w = self.stencil(slice(xiis,xiie),xii)
				for dn in dns:
					bcw = self.bcleft.jacobian_left(w[dn])
					if bcw.ndim==1:
						self.dmatrix[dn].lines[ii] = Line(bcw, ii)
					else:
						self.dmatrix[dn].lines[ii] = LineMulti(bcw, ii)

		if rightbc:
			#Right boundary conditions
			for ii in range(Nx-bcend,Nx):
				xii = ii+bco
				xiis = xii-ox; xiie = self.Nx_bc
	
				#Apply boundary conditions to stencil based on x with boundary
				#nodes included, but centered at ii+bcoffset
				w = self.stencil(slice(xiis,xiie),xii)
				for dn in dns:
					bcw = self.bcright.jacobian_right(w[dn])
					if bcw.ndim==1:
						self.dmatrix[dn].lines[ii] = Line(bcw, offsetr=Nx-ii)
					else:
						self.dmatrix[dn].lines[ii] = LineMulti(bcw, offsetr=Nx-ii)
		
		#If we only have the final BC then this would speed it up
		#for dn in dns:
		#	self.dmatrix[dn].line_range = (Nx-1,Nx)
		
	def diff(self, x, n=1):
		"Apply the Jacobian matrix J_n to vector x"
		if n==0: return x*0
		return self.dmatrix[n].matvec(x)

	def rdiff(self, x, n=1):
		"Apply the adjoint of the Jacobian matrix (J_n)^A to vector x"
		if n==0: return x*0
		return self.dmatrix[n].rmatvec(x)

