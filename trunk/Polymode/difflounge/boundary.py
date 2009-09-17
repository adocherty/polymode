# _*_ coding=utf-8 _*_

'''
Changes:
 	* Implemented n element dtn functions
	* Changed BC to have apply_left & apply_right members
	* Left BC .. with Dirichlet/Neumann BC
Todo:
	* Include bcorder for different orders of boundary conditions
'''
from __future__ import division
import logging

from numpy import *
from .finitedifference import fd_calc_weights

aa_=asarray

## Boundary condtions!!
class BoundaryCondition:
	bcorder=0
	row_extent=1
	def apply_left(self, w, row):
		return w
	def apply_right(self, w, row):
		return w

	#Update stencil
	def set_stencil(self, x):
		self.bc_location_left = x[0]
		self.bc_location_right = x[-1]

	#BC with adjustable paramters
	def set_condition(self, c):
		pass

	def __str__(self):
		return self.__doc__

BC_None=BoundaryCondition

## Dirichlet BC	
class BC_Dirichlet(BoundaryCondition):
	"Dirichlet Boundary Conditions"
	def __init__(self, X=None, dtype=float_):
		self.dtype = dtype
		self.bcorder = 1
		if iterable(X): self.set_stencil(X)

	def apply_left(self, w, row=0):
		return w[1:]

	def apply_right(self, w, row=0):
		return w[:-1]

	def jacobian_left(self, w, row=0):
		return 0*w[1:]
		
	def jacobian_right(self, w, row=0):
		return 0*w[:-1]


#Mixed/Neumann Boundary Condition
class BC_Mixed(BoundaryCondition):
	"""
	General Mixed Boundary Conditions
	dy/dx + D y = 0
	"""
	def __init__(self, X=None, bandwidth=None, xbc=0, condition=0, dtype=float_):
		self.bandwidth = bandwidth
		self.dtype = dtype
		self.xbc = xbc
		self.set_condition(condition)

		#Set the bc stencil
		if X is not None:
			self.set_stencil(X)
		
		#used by the FD routine to know how many boundary points we have eliminated
		#Currently only works for 1 - standard 2nd order BCs
		self.bcorder = 1

		#Offset for the finite difference stencil
		self.offset = bandwidth//2

	def __str__(self):
		outstr = self.__doc__+", bw: %d, xbc: %.2g" % (self.bandwidth, self.xbc)
		return outstr

	def set_condition(self,dtn,cf=1,cm=1):
		self.dtn=asarray(dtn)
		self.chain_factor = cf
		self.cm = cm

	#Overloadable function to calculate bc stencil
	def set_stencil(self, x):
		xbcl = self.bc_location_left = self.xbc*(x[0]-x[1])+x[1]
		xbcr = self.bc_location_right = self.xbc*(x[-1]-x[-2])+x[-2]
		
		bw = self.bandwidth
		self.bc_stencil_left = fd_calc_weights(x[:bw],xbcl,1,dtype=self.dtype)
		self.bc_stencil_right = fd_calc_weights(x[-bw:],xbcr,1,dtype=self.dtype)

	def apply_left(self, w, row=0):
		bcs = self.bc_stencil_left
		bcs = bcs.reshape(bcs.shape+(1,)*self.dtn.ndim)
		bcshape = (max(self.bandwidth,w.shape[0]),) + shape(self.dtn)

		rw = zeros(bcshape, dtype=self.dtype)
		rw.T[...,:w.shape[-1]] = w*self.cm
		rw[:self.bandwidth] +=  w[0]*self.dtn[newaxis]*(bcs[0]*bcs[1][0]-bcs[1]*bcs[0][0]) \
									/(bcs[1][0]-self.dtn[newaxis]*bcs[0][0])/bcs[1][0]
		rw[:self.bandwidth] +=  -w[0]*self.cm*bcs[1]/bcs[1][0]
		return rw[1:]

	def apply_right(self, w, row=0):
		bcs = self.bc_stencil_right
		bcs = bcs.reshape(bcs.shape+(1,)*self.dtn.ndim)
		bcshape = (max(self.bandwidth,w.shape[0]),) + shape(self.dtn)
		
		rw = zeros(bcshape, dtype=self.dtype)
		rw.T[...,-w.shape[0]:] = w*self.cm
		rw[-self.bandwidth:] += w[-1]*self.dtn[newaxis]*(bcs[0]*bcs[1][-1]-bcs[1]*bcs[0][-1]) \
							/(bcs[1][-1]-self.dtn[newaxis]*bcs[0][-1])/bcs[1][-1]
		rw[-self.bandwidth:] += -w[-1]*self.cm*bcs[1]/bcs[1][-1]
		return rw[:-1]

	def jacobian_left(self, w, row=0):
		bcs = self.bc_stencil_left
		bcs = bcs.reshape(bcs.shape+(1,)*self.dtn.ndim)
		bcshape = (max(self.bandwidth,w.shape[0]),) + shape(self.dtn)

		jac = zeros(bcshape, dtype=self.dtype)
		jac[:self.bandwidth] = w[0]*(bcs[1][0]*bcs[0]-bcs[0][0]*bcs[1]) \
								/(bcs[1][0]-self.dtn[newaxis]*bcs[0][0])**2
		return jac[1:]*self.chain_factor
		
	def jacobian_right(self, w, row=0):
		bcs = self.bc_stencil_right
		bcs = bcs.reshape(bcs.shape+(1,)*self.dtn.ndim)
		bcshape = (max(self.bandwidth,w.shape[0]),) + shape(self.dtn)

		jac = zeros(bcshape, dtype=self.dtype)
		jac[-self.bandwidth:] = w[-1]*(bcs[1][-1]*bcs[0]-bcs[0][-1]*bcs[1]) \
								/(bcs[1][-1]-self.dtn[newaxis]*bcs[0][-1])**2
		return jac[:-1]*self.chain_factor


#Mixed/Neumann Boundary Condition
class BC_OnePoint(BoundaryCondition):
	"""
	General Mixed Boundary Conditions
	dy/dx + D y = 0
	involving only the first & last points of the function
	"""
	def __init__(self, X=None, condition=0, dtype=float_):
		self.dtn = condition
		self.dtype = dtype
		self.set_condition(condition)
		
		#Set the bc stencil
		if X is not None:
			self.set_stencil(X)
		
		#used by the FD routine to know how many boundary points we have eliminated
		#Currently only works for 1 - standard 2nd order BCs
		self.bcorder = 1
		
	def set_stencil(self, x):
		self.bc_location_left = 0.5*(x[1]+x[0])
		self.bc_location_right = 0.5*(x[-1]+x[-2])
		self.h_left = x[1]-x[0]
		self.h_right = x[-1]-x[-2]
	
	def set_condition(self,dtn,cf=1):
		self.dtn=asarray(dtn)
		self.chain_factor = cf

	def apply_left(self, w, row=0):
		rw = zeros(w.shape + shape(self.dtn), dtype=self.dtype)
		rw.T[:] = w
		rw[1] += (2-self.dtn*self.h_left)/(2+self.dtn*self.h_left)*w[0]
		return rw[1:]

	def apply_right(self, w, row=0):
		rw = zeros(w.shape + shape(self.dtn), dtype=self.dtype)
		rw.T[:] = w
		rw[-2] += (2+self.dtn*self.h_right)/(2-self.dtn*self.h_right)*w[-1]
		return rw[:-1]

	def jacobian_left(self, w, row=0):
		rw = zeros(w.shape + shape(self.dtn), dtype=self.dtype)
		rw[1] = -4*self.h_right/(2+self.dtn*self.h_right)**2*w[0]
		return rw[1:]*self.chain_factor
		
	def jacobian_right(self, w, row=0):
		rw = zeros(w.shape + shape(self.dtn), dtype=self.dtype)
		rw[-2] = 4*self.h_right/(2-self.dtn*self.h_right)**2*w[-1]
		return rw[:-1]*self.chain_factor


class BC_Neumann(BC_Mixed):
	"Neuman Boundary Conditions"
	def __init__(self, X=None, bandwidth=3, xbc=0, dtype=float_):
		BC_Mixed.__init__(self, X=X, bandwidth=bandwidth, xbc=xbc, condition=0, dtype=dtype)

	def set_condition(self,dtn,cf=1):
		self.dtn=asarray(dtn)*0
		self.chain_factor=cf


class BC_Switch(BC_Neumann):
	"Switched Dirichlet/Neuman Boundary Conditions"
	def __init__(self, X=None, bandwidth=3, xbc=1, mask=0, dtype=complex):
		self.mask = array(mask)
		BC_Mixed.__init__(self, X=X, bandwidth=bandwidth, condition=0, xbc=xbc, dtype=dtype)
	
	def set_mask(self,mask):
		self.mask = asarray(mask)
		
	def apply_left(self, w, row=0):
		bcs = self.bc_stencil_left
		bcs = bcs.reshape(bcs.shape+(1,)*self.mask.ndim)
		rw = zeros((max(self.bandwidth,w.shape[0]),) + shape(self.mask), dtype=self.dtype)
	
		rw.T[...,:w.shape[0]] = w
		rw[:self.bandwidth] +=  self.mask[newaxis]*w[0]*(-bcs[1])/(bcs[1][0])
		return rw[1:]

	def apply_right(self, w, row=0):
		bcs = self.bc_stencil_right
		bcs = bcs.reshape(bcs.shape+(1,)*self.mask.ndim)
		rw = zeros((max(self.bandwidth,w.shape[0]),) + shape(self.mask), dtype=self.dtype)
	
		rw.T[...,-wshape[0]:] = w
		rw[...,-self.bandwidth:] += self.mask[newaxis]*w[-1]*(-bcs[1])/(bcs[1][-1])
		return rw[:-1]


