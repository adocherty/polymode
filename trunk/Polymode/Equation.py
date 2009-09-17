# _*_ coding=utf-8 _*_
#---------------------------------------------------------------------------------
#Copyright © 2009 Andrew Docherty
#
#This program is part of Polymode.
#Polymode is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#---------------------------------------------------------------------------------
"""
The action construction for the Vector Wave Equation

Todo:
 - Extend to use finite differences in the azimuthal direction
 - Add NRBC at the inner radius to enable reduced sector computations
 - Add reflection symmetry
"""

from __future__ import division
import logging

#Numpy and Scipy imports
from numpy import *
from numpy.lib.scimath import sqrt

#Helper library imports
from .mathlink.blockarray import BlockArray
from .mathlink import timer, fftp, hankel1, hankel1p, hankel2, hankel2p
from .mathlink import hankel1_ratio, besselj_ratio

#ABC Solver imports
from . import Modes

class Equation(object):
	pass
	
##Action of SWE equation on vector DM in Fourier domain
class VectorWaveEquation(Equation):
	"""
	VWE Magnetic equation class:
	diff: Difference class for finite differences
	Calculate the action of the VWE on a vector, defined as
	L^p H^p + G^{p,q=1} H^+ + G^{p,q=-1} H^-
	for x^p
	
	Note: x must be shaped correctly, as x[Np, Nr, bw, Naz, Naz] with Np=2
	"""
	r_axis = 0
	phi_axis = 1
	p_axis = 2
	
	def __init__(self, diff, fast_convolve=1, modetype='outward', dtype=complex_):
		"""
		Create a VWE Magnetic equation object
		
		Parameters:
		  -  `dtype`: Numpy data type
		  -  `diff`: Difference class to used internally
		  -  `modetype`: 'outwards', 'inwards' or 'bound'
		"""
		self.dtype = dtype
		self.diff = diff
		self.modetype = modetype
		
		#the number of polarization components
		self.pmax = 2
		self.pv = array([1,-1])

		#Paramters
		self.fftnum = 0
		self.fast_convolve = fast_convolve

	def __str__(self):
		infostr = "Vector Wave Equation with "
		infostr += u"m0=%d, \u03bb=%.4g\n".encode('utf-8') % (self.m0, 2*pi/self.k0)
		infostr += self.diff.__str__()
		return infostr

	def __getstate__(self):
		"Pickle all needed data, ignore cached data"
		state = self.__dict__.copy()
		ignore_list = ["ni2", "drlogn", "dazlogn"]
		for ignore in ignore_list:
			if ignore in state: del state[ignore]
		return state
	
	def __setstate__(self,state):
		"Restore pickled data"
		self.__dict__.update(state)

	def setup(self,shape,wg,m0,wl):
		self.m0 = m0
		self.wavelength = wl
		self.k0 = 2*pi/wl
		self.shape = shape + (self.pmax,)
		self.fd_shape = self.diff.bandwidth,shape[1],self.pmax

		self.wg = wg
		self.exterior_index = wg.exterior_index(wl)
		self.interior_index = wg.interior_index(wl)
		
		#This should be more adjustable
		#Change it when updating to new azimuthal difference class
		if self.fast_convolve or shape[1]==1:
			self.wgshape = shape[0], shape[1]
		else:
			self.wgshape = shape[0], 2*shape[1]	
		
		#Coordinate information
		wgcoord = self.wg.get_coord(self.wgshape, border=1)
		self.coord = wg.get_coord(shape, border=1)
		self.msm0 = (self.coord.ms[:,newaxis]+self.m0)

		#Construct the finite differences (BCs are updated later)
		self.diff.generate()

		#Index and derivatives
		cinx = self.wg.calculate(wgcoord, self.wavelength)
		
		#self.ni2 = fftp.ifftshift(cinx.ni2, axes=[1])[...,newaxis]
		#self.drlogn = fftp.ifftshift(cinx.drlogn, axes=[1])[...,newaxis]
		#self.dazlogn = fftp.ifftshift(cinx.dazlogn, axes=[1])[...,newaxis]
		self.ni2 = cinx.ni2[...,newaxis]
		self.drlogn = cinx.drlogn[...,newaxis]
		self.dazlogn = cinx.dazlogn[...,newaxis]
		
	#Linear / Circular convolution
	def fftmultiply(self, x,fy, axis=1):
		if self.fast_convolve:
			y = fftp.ifft(fy, axis=axis, overwrite_x=0)
			y *= x
			fftp.fft(y, axis=axis, overwrite_x=1)
			return y
		else:
			Naz = self.shape[axis]

			fyb = zeros((fy.shape[0], x.shape[1], fy.shape[2]), fy.dtype)	
			xy = zeros(fy.shape, fy.dtype)

			fyb[:,:(Naz+1)//2] = fy[:,:(Naz+1)//2]
			fyb[:,-(Naz//2):] = fy[:,-(Naz//2):]
			xyb = fftp.fft(x*fftp.ifft(fyb, axis=1), axis=1)
			
			#Resize
			xy[:,:(Naz+1)//2] = xyb[:,:(Naz+1)//2]
			xy[:,-(Naz//2):] = xyb[:,-(Naz//2):]
			return xy

	def calc_rightbc_dtn(self, ev, rdtn=None):
		"Calculate Dirichlet-to-Neumann coefficients from eigenvalue"
		g = Modes.branchsqrt((self.exterior_index*self.k0)**2-ev)
		if rdtn is None: rdtn = self.diff.bcright.bc_location_right
				
		#For g~0 set Neumann BC, otherwise we get NaN's
		if isnan(g) or abs(ev)>1e4 or abs(g)<1e-10:
			bcdtn = zeros((self.coord.ms.shape[0], self.pv.shape[0]))
		else:
			bcdtn = hankel1_ratio(self.msm0+self.pv, g*rdtn)/rdtn
		return bcdtn

	def calc_leftbc_dtn(self, ev, rdtn=None):
		"Calculate Dirichlet-to-Neumann coefficients from eigenvalue"
		g = Modes.branchsqrt((self.interior_index*self.k0)**2-ev)
		if rdtn is None: rdtn = self.diff.bcleft.bc_location_left
		
		bcdtn = besselj_ratio(self.msm0+self.pv, g*rdtn)/rdtn
		return bcdtn

	def set_lambda(self, ev, onlyev=0, onlyconst=0):
		"Set the eigenvalue for the boundary conditions"
		if self.coord.rmin==0:
			self.diff.bcleft.set_mask((self.msm0+self.pv)==0)
			if onlyconst:
				self.diff.bcright.set_condition(0, cm=1)
			elif onlyev:
				self.diff.bcright.set_condition(self.calc_rightbc_dtn(ev), cm=0)
			else:
				self.diff.bcleft.set_condition(self.calc_rightbc_dtn(ev))
		else:
			if onlyconst:
				self.diff.bcleft.set_condition(0, cm=1)
				self.diff.bcright.set_condition(0, cm=1)
			elif onlyev:
				self.diff.bcleft.set_condition(self.calc_leftbc_dtn(ev), cm=0)
				self.diff.bcright.set_condition(self.calc_rightbc_dtn(ev), cm=0)
			else:
				self.diff.bcleft.set_condition(self.calc_leftbc_dtn(ev))
				self.diff.bcright.set_condition(self.calc_rightbc_dtn(ev))
		self.diff.generate(leftbc=1, rightbc=1)

	def Gx(self, x, r, drlogn, dazlogn, D, p=1, q=1):
		Gl = p*drlogn + 1j*dazlogn/r
		G = -0.5*self.fftmultiply(Gl,q*D.diff1(x)+(self.msm0+q)*D.diff0(x)/r )
		return G

	def Gxr(self, x, r, drlogn, dazlogn, D, p=1, q=1):
		gx = conj(self.fftmultiply(conj(p*drlogn + 1j*dazlogn/r),conj(x)))
		G = -0.5*(q*D.rdiff1(gx) + (self.msm0+q)*D.rdiff0(gx)/r)
		return G

	def Lx(self, x, r, ni2, D, p=1):
		Lx = D.diff2(x) + D.diff1(x)/r + D.diff0(-(self.msm0+p)**2/r**2*x \
			+ self.k0**2*self.fftmultiply(ni2, x))
		return Lx

	def Lxr(self, x, r, ni2, D, p=1):
		Lx = D.rdiff2(x) + D.rdiff1(x/r) + D.rdiff0(-(self.msm0+p)**2/r**2*x \
			+ self.k0**2*conj(self.fftmultiply(conj(ni2), conj(x))))
		return Lx

	def matvec(self, x):
		"Apply the VWE to the vector x"
		xshape = x.shape
		if x.ndim<>3: x = x.reshape(self.shape)
		r = self.coord.rv[:, newaxis, newaxis]
		
		#The SWE L and the vector contribution G and cross-contribution Gx
		L = self.Lx(x, r, self.ni2, self.diff, self.pv)
		G = self.Gx(x, r, self.drlogn, self.dazlogn, self.diff, self.pv, self.pv)
		Gx = self.Gx(x, r, self.drlogn, self.dazlogn, self.diff, -self.pv, self.pv)

		return (L+G+Gx[...,::-1]).reshape(xshape)

	def rmatvec(self, x):
		"Apply the adjoint of the VWE to the vector x"
		xshape=x.shape
		if x.ndim<>3: x = x.reshape(self.shape)
		x = conj(x)
		r = self.coord.rv[:, newaxis, newaxis]

		#The SWE L and the vector contribution G and cross-contribution Gx
		L = self.Lxr(x, r, self.ni2, self.diff, self.pv)
		G = self.Gxr(x, r, self.drlogn, self.dazlogn, self.diff, self.pv, self.pv)
		Gx = self.Gxr(x[...,::-1], r, self.drlogn, self.dazlogn, self.diff, -self.pv, self.pv)

		return conj(L+G+Gx).reshape(xshape)

	def construct(self, x, row=0):
		"Apply the VWE to construct the matrix"
		r = self.coord.rv[row:row+1, newaxis, newaxis]
		drln = self.drlogn[row:row+1]
		daln = self.dazlogn[row:row+1]
		
		#The SWE L and the vector contribution G and cross-contribution Gx
		D = self.diff.construct(row)
		L = self.Lx(x, r, self.ni2[row:row+1], D, self.pv)
		G = self.Gx(x, r, drln, daln, D, self.pv, self.pv)
		Gx = self.Gx(x, r, drln, daln, D, -self.pv, self.pv)

		return L+G+Gx[...,::-1]

	def __call__(self, x):
		return self.matvec(x)

class VectorWaveJacobian(VectorWaveEquation):
	def set_lambda(self, ev):
		"Set the eigenvalue for the boundary conditions"
		#Caclulate the chain rule factor ∂D/∂λ
		mm = self.coord.ms + self.m0

		g = Modes.branchsqrt((self.exterior_index*self.k0)**2-ev)
		rc = self.diff.bcright.bc_location_right
		bcfacp_right = 0.5*rc*array([hankel1_ratio(mm+1,g*rc)**2/(g*rc)**2 + (1-(mm+1)**2/(g*rc)**2),
					hankel1_ratio(mm-1,g*rc)**2/(g*rc)**2 + (1-(mm-1)**2/(g*rc)**2)]).T

		g = sqrt(self.interior_index**2*self.k0**2 - ev)
		rc = self.diff.bcleft.bc_location_left
		bcfacp_left = 0.5*rc*array([besselj_ratio(mm+1,g*rc)**2/(g*rc)**2 + (1-(mm+1)**2/(g*rc)**2),
					besselj_ratio(mm-1,g*rc)**2/(g*rc)**2 + (1-(mm-1)**2/(g*rc)**2)]).T

		#Now calculate the standard Jacobian BCs
		self.diff.set_left_bc(self.calc_leftbc_dtn(ev), bcfacp_left)
		self.diff.set_right_bc(self.calc_rightbc_dtn(ev), bcfacp_right)
		self.diff.generate(leftbc=1, rightbc=1)

