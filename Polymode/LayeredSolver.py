# _*_ coding=utf-8 _*_
#
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

from __future__ import division
import logging

import pylab as pl
from numpy import *
from scipy import sqrt,optimize

from .Solver import *
from .Modes import Mode, branchsqrt
from .mathlink import hankel1, hankel2, hankel1p, hankel2p, jv, jvp
from .mathlink import coordinates,timer

def hankel1p(m, x):
	return 0.5*(hankel1(m-1,x)-hankel1(m+1,x))
def hankel2p(m, x):
	return 0.5*(hankel2(m-1,x)-hankel2(m+1,x))

def convert_condon_to_dbf(nc, gammac, k0, V=1.0, tol=1e-8):
	"Convert Condon chiral paramters to Drude-Born-Federov"
	Gdbf=1
	change=inf; ii=0
	while (change>tol) or (ii>5):
		n_dbf =nc*Gdbf
		gamma_dbf = gammac*V/(k0*n_dbf**2)*Gdbf
		k = k0*n_dbf
		Gdbf = 1-k**2*gamma_dbf**2

		change = abs(gamma_dbf-gammac*V/(k0*n_dbf**2)*Gdbf).max()

	return n_dbf, gamma_dbf

def find_sign_change(x, n=1, axis=-1):
	"""
	Find sign changes of the data along axis
	"""
	#Find signs
	schange = sign(x)
	
	#Iterate over all other axes
	scindex = transpose(diff(schange, n=n, axis=axis).nonzero())
	return scindex

class Layer(object):
	'''
	The base class for a layer in the transfer matrix approach
	'''
	shape = (4,4)

	def __init__(self, r1, r2, n, gammac=0, zorder=0):
		self.r1,self.r2,self.zorder = r1,r2,zorder
		self.n,self.gammac = n,gammac

	def __cmp__(self, other):
		#Compare two layers
		if hasattr(other,'r1'):
			return cmp(self.r1,other.r1)
		else:
			return cmp(self.r1,other)
	
	def __str__(self):
		return "%s: r1=%g r2=%g n=%g z=%g" % (type(self).__name__, self.r1, self.r2, self.n, self.zorder)
	def __repr__(self):
		return "%s(%g,%g,%g,%g,%g)" % (type(self).__name__, self.r1, self.r2, self.n, self. gammac, self.zorder)

	def copy(self):
		other = self.__class__(self.r1,self.r2,self.n,self.gammac,self.zorder)
		other.setup(self.m, self.k0)
		return other
				
	def setup(self, m, k0, coeffs=None):
		self.coeffs = coeffs
		self.m = m
		self.k0 = k0
		self.al = 1j*self.n
		self.ar = -1/self.al
	
		k = k0*self.n
		self.kp = k/(1-self.gammac*k)
		self.km = k/(1+self.gammac*k)
		self.Kpm = array([[-self.km,0],[0,self.kp]])
		self.V = array([[self.al,1],[1,self.ar]])

	def root(self,x):
		return branchsqrt(x)

	def precalc(self, beta):
		Vm = self.km**2 - beta**2
		Vp = self.kp**2 - beta**2
		xim, xip = self.root(Vm), self.root(Vp)
		S = array([[Vp/Vm, 0],[0,1]])
		
		self.precalc_data = Vm,Vp,xim,xip,S
		return self.precalc_data

	def multiplicative_factor(self):
		Vm,Vp,xim,xip,S = self.precalc_data
		return Vp
		
	def calculate(self, beta, r, precalc=1):
		kp,km = self.kp,self.km
		
		#Precalculate beta dependant stuff to speed things up
		if precalc:
			Vm,Vp,xim,xip,S = self.precalc(beta)
		else:
			Vm,Vp,xim,xip,S = self.precalc_data
	
		#Calculate transfer matrix
		Q,dQ = self.Q(xim,xip,r)
		T = zeros(self.shape, dtype=complex_)
		T[:2, :Q.shape[1]] = dot(self.V,Q)*Vp
		T[2:, :Q.shape[1]] = -dot(dot(self.V,S),(beta*self.m*Q/r + dot(self.Kpm,dQ)))
		return T

	def field(self, beta, r, he=True, precalc=True):
		if precalc:
			Vm,Vp,xim,xip,S = self.precalc(beta)
		else:
			Vm,Vp,xim,xip,S = self.precalc_data

		Q,dQ = self.Q(xim,xip,r)

		#Multiply by coefficients
		av = self.coeffs[:Q.shape[1]]
		qz = self.multiplicative_factor()*dot(Q,av)
		
		#Deal with different cases as r->0
		if abs(r)>1e-10:
			qa = -dot(dot(S,beta*self.m*Q/r + dot(self.Kpm,dQ)),av)
			qr = 1j*dot(dot(S,dot(self.Kpm,Q)*self.m/r + beta*dQ),av)
		elif self.m==0:
			qa = -dot(dot(S,dot(self.Kpm,dQ)),av)
			qr = 1j*dot(dot(S,beta*dQ),av)
		else:
			qa = -dot(dot(S,self.m*beta*dQ + dot(self.Kpm,dQ)),av)
			qr = 1j*dot(dot(S,self.m*dot(self.Kpm,dQ) + beta*dQ),av)
		
		#Convert to magnetic & electric fields
		if he:
			F = array([dot(self.V,qz),dot(self.V,qa),dot(self.V,qr)]).T
		else:
			F = array([qz,qa,qr]).T

		return F
	
class HLayer(Layer):
	def Q(self, xim, xip, r):
		m = self.m
		Z = zeros_like(r)
		#Define the base Q
		Q = array([[hankel1(m,xim*r), Z, hankel2(m,xim*r), Z],
					[Z, hankel1(m,xip*r), Z, hankel2(m,xip*r)]])
		#Derivative of Q
		dQ = array([[xim*hankel1p(m,xim*r), Z, xim*hankel2p(m,xim*r), Z],
					[Z, xip*hankel1p(m,xip*r), Z, xip*hankel2p(m,xip*r)]])
		return Q,dQ

class HLayerExterior(HLayer):
	def Q(self, xim, xip, r):
		m = self.m
		Z = zeros_like(r)
		#Define the base Q
		Q = array([[hankel1(m,xim*r), Z], [Z, hankel1(m,xip*r)]])
		#Derivative of Q
		dQ = array([[xim*hankel1p(m,xim*r), Z], [Z, xip*hankel1p(m,xip*r)]])
		return Q,dQ

class KLayerExterior(Layer):
	def root(self,x):
		return branchsqrt(-x)
	def Q(self, xim, xip, r):
		m = self.m
		Z = zeros_like(r)
		#Define the base Q
		Q = array([[kv(m,xim*r), Z], [Z, kv(m,xip*r)]])
		#Derivative of Q
		dQ = array([[xim*kvp(m,xim*r), Z], [Z, xip*kvp(m,xip*r)]])
		return Q,dQ

class JLayerInterior(Layer):
	def Q(self, xim, xip, r):
		m = self.m
		Z = zeros_like(r)
		#Define the base Q
		Q = array([[jv(m,xim*r), Z], [Z, jv(m,xip*r)]])
		#Derivative of Q
		dQ = array([[xim*jvp(m,xim*r), Z], [Z, xip*jvp(m,xip*r)]])
		return Q,dQ

InteriorLayer=JLayerInterior
MidLayer=HLayer
ExteriorLayer=HLayerExterior

class LayeredMode(Mode):
	def normalize(self, by='power', wg=None, coord=None):
		"""
		Normalize the fields so the electric and magnetic vectors have the correct
		correspondance, including the intrinsic impedence of free space:
		H = Htrue, E = (e0/mu0)^1/2 Etrue
		"""
		pass

	def magnetic_transverse_field(self, fourier=False, cartesian=None, coord=None):
		'''
		The transverse magnetic field, calculated from the internal H⁺,H⁻"
		cartesian=False returns h_t=(h_r,h_ϕ)
		cartesian=True returns h_t=(h_x,h_y)
		'''
		H,E=self._construct_fields_(coord, he=True)
		return H[:2]

	def magnetic_field(self, fourier=False, cartesian=None, coord=None):
		"""
		The three component magnetic field (h_r,h_ϕ,h_z) or (hx,hy,hz)
		"""
		H,E=self._construct_fields_(coord, he=True)
		return H

	def electric_transverse_field(self, fourier=False, cartesian=None, coord=None, wg=None):
		'''
		The transverse electric field, calculated from the internal E⁺,E⁻"
		cartesian=False returns e_t=(e_r,e_ϕ)
		cartesian=True returns e_t=(e_x,e_y)
		
		if calculated_electric_field is true then calculate from the the magnetic field
		'''
		if coord is None: #coord=self.coord
			raise NotImplementedError, "Need a coord for LayeredMode"
		
		H,E=self._construct_fields_(coord, he=True)
		return E[:2]

	def electric_field(self, wg=None, fourier=False, cartesian=None, coord=None):
		"""
		The three component electric field (e_r,e_ϕ,e_z)
		if calculated_electric_field is true then calculate from the the magnetic field
		"""
		H,E=self._construct_fields_(coord, he=True)
		return E

	def _construct_fields_(self, coord, he=True):
		"Sample fields on a coord grid"
		beta = self.beta

		#Get r, phi pairs for all points
		rm, phim = coord.polar2d()

		#Sort by rm - we calculate the fields for each r only once
		sri = argsort(rm.flat)
		rmflat = rm.flat[sri]
	
		#Calculate the fields
		F = zeros((2,3,len(rmflat)), dtype=complex_)
		fshape = (2,3)+rm.shape
		
		#Iterate over each layer
		Nlayer = len(self.layers)
		rmin = 0
		for ii in xrange(0, Nlayer):
			layer = self.layers[ii]

			#Select points inside layer
			lstart = rmflat.searchsorted(layer.r1)
			lend = rmflat.searchsorted(layer.r2)

			#all rm within this layer
			number_different=0; last_r = None
			for kk in xrange(lstart,lend):
				sii = sri[kk]; r = rm.flat[sii]; phi = phim.flat[sii]
				
				#Check if r change changed
				if r!=last_r:
					Fcurrent = layer.field(beta, r, he=he)
				
				F[...,kk] = Fcurrent*exp(1j*self.m0*phi)
		return F.reshape(fshape)

class LayeredSolver(Solve):
	def setup(self, Nscan=1e3, tol=1e-6):
		self.Nscan = Nscan
		self.tolerance = tol

	def setup_layers(self):
		#Construct layers
		self.layers = self.wg.calculate_radial_layers(self.wl)
		self.Nlayer = len(self.layers)

		#Setup other parameters -- a bit of a hack!
		for layer in self.layers:
			layer.setup(self.m0, self.k0)

		self.C = zeros((4,4), dtype=complex_)
		self.J = zeros((4,4), dtype=complex_)

	def krange(self):
		"Returns range of k scanned over, real and imag"
		#Set krange automatically or manually
		if self.bracket is None:
			#Set automatic nrange unless specific ns given
			kmin = min([min(l.kp,l.km) for l in self.layers])
			kmax = max([max(l.kp,l.km) for l in self.layers])
			krange = array([kmin,kmax,0,0] if self.nefflist is None else [0,inf,0,0])

		else:
			#Set manual nrange
			krange = append(real(self.bracket), imag(self.bracket))*self.k0
		return krange

	def nmax(self):
		nplus = max([l.kp/self.k0 for l in self.layers])
		nminus = max([l.km/self.k0 for l in self.layers])
		return nplus,nminus
		
	def set_lambda(self, beta):
		self.get_matrix(beta)
		
	@timer.time_function(prefix="T")
	def get_matrix(self, beta):
		"Return the matrix for the modal eigenvalue problem"

		#Initial layer coefficient matrix
		To = self.layers[0].calculate(beta,self.layers[0].r2)[:,:2]

		#Final layer coefficient matrix
		Tn = self.layers[-1].calculate(beta,self.layers[-1].r1)[:,:2]
		
		#Intermediate layers
		for ii in range(1, self.Nlayer-1):
			layer = self.layers[ii]

			T1 = layer.calculate(beta,layer.r1)
			T2 = layer.calculate(beta,layer.r2, precalc=0)
			To = dot(T2,linalg.solve(T1,To))

		#V factor - better way to implement this?
		Tn /= self.layers[-1].multiplicative_factor()
		To /= self.layers[0].multiplicative_factor()

		#Assert conditions on first an last layer
		self.C[:, :2] = To[:,:2]
		self.C[:, 2:] = -Tn[:,:2]
		return self.C

	def get_jacobian(self, beta):
		pass

	def det(self, betas):
		#betas = squeeze(betas)
		if iterable(betas):
			ans = zeros(shape(betas), complex)
			for inx in ndindex(*shape(betas)):
				ans[inx] = self.det(betas[inx])
		else:
			ans = linalg.det(self.get_matrix(betas))
		return ans

	def det_ri(self, bri):
		"Interface for complex root finder"
		beta = bri[0]+bri[1]*1j
		det = self.det(beta)
		return [det.real, det.imag]

	@timer.time_function()
	def calculate(self, number=inf):
		"""
		Scan for modes
		"""
		roots = []
		krange = self.krange()
		k0 = self.k0
		scan_complex = abs(krange[-1]-krange[-2])>0
		logging.info( "Searching range of neffs: %.5g -> %.5g" % (krange[0]/k0, krange[1]/k0) )

		#Create layers from waveguide
		self.setup_layers()
		modecoord = coordinates.PolarCoord(rrange=(0,self.layers[-1].r1), arange=(-pi,pi))

		#Iterate on list OR search over krange
		if self.nefflist is not None:
			kscan = atleast_1d(self.nefflist)*k0
			possible_zc = transpose(kscan.nonzero())
		
		elif self.modelist is not None:
			nefflist = [m.neff for m in self.modelist]
			kscan = atleast_1d(nefflist)*k0
			possible_zc = transpose(kscan.nonzero())
		
		else:
			#Define course mesh for scanning
			dkr = (krange[1]-krange[0])/self.Nscan
			kr = arange(real(krange[0])+dkr,real(krange[1])-dkr,dkr)

			#Scan over complex neff if leaky
			if scan_complex:
				dki = abs(krange[3]-krange[2])/5
				ki = arange(krange[2], krange[3], dki)
				kscan = kr + ki[:,newaxis]*1j
			else:
				kscan = kr

			logging.debug( "Scanning %d points" % len(kscan.flat) )

			#Calculate determinant of transfer matrix over range
			tix = self.det(kscan)
	
			#Detect local minima/maxima
			possible_zc = find_sign_change(diff(absolute(tix)))[::-1]

		#Iterate over all possible zeros and seek closest zero
		logging.info( "Searching for %d possible modes" % len(possible_zc) )
		kk = 0
		for inx in possible_zc:
			inx = tuple(inx)
			#Try and locate closest root, if we fail go to the next in the list
			try:
				bri  = optimize.fsolve(self.det_ri, [kscan[inx].real,kscan[inx].imag], warning=False, xtol=1e-12)
				root = complex(*bri)
			except linalg.LinAlgError:
				continue

			#ignore those outside the range
			if real(root)<krange[0] or real(root)>krange[1]:
				logging.debug( "Rejecting out of range solution" )
				continue

			#Ignore those with a large residue
			res = abs(self.det(root))
			if res>self.tolerance:
				logging.debug( "Rejecting inaccurate solution" )
				continue

			#Finally add the root if not already in the list
			if (len(roots)==0 or min(abs(root-array(roots)))>self.tolerance):
				mode = LayeredMode(m0=self.m0, wl=self.wl, coord=modecoord, symmetry=1, evalue=root**2)
				mode.layers = self.calculate_mode_layers(root)
				mode.residue = res
				self.modes += [mode]

				logging.info( "Mode #%d: neff=%s, res: %.2e\n" % (kk, root/k0, res) )
				roots.append(root)
				kk+=1
				self.numbercalculated+=1
					
			if kk>=number or (self.numbercalculated>=self.totalnumber):
				break
	
		#Sort modes in finalization method
		self.modes.sort()
		return self.modes

	def calculate_mode_layers(self, beta):
		"Create new layers list with coefficients"
		#Find eigenvector(s)
		A = self.get_matrix(beta)
		w,v = linalg.eig(A)
		
		mode_inx = abs(w).argmin()
		#if len(mode_inx)>1:
		#	print "Found %d degenerate modes" % len(mode_inx)
		assert abs(w[mode_inx])<1e-4, "Error in field calculation: mode not accurate"

		#Field coefficients
		mode_coeffs = v[:,mode_inx] #/absmax(v[:,mode_inx])
		a0 = zeros_like(mode_coeffs)
		an = zeros_like(mode_coeffs)

		a0[:2] = mode_coeffs[:2]
		an[:2] = mode_coeffs[2:]

		#Create new layer list with coefficients
		last = self.layers[0].copy()
		last.coeffs = a0
		mode_layers = [last]

		for ii in range(1, self.Nlayer-1):
			layer = self.layers[ii].copy()
			T1 = last.calculate(beta, layer.r1)
			T2 = layer.calculate(beta, layer.r1)

			layer.coeffs = dot(linalg.inv(T2),dot(T1,last.coeffs))
			mode_layers.append(layer)

		#Final layer:
		layer = self.layers[-1].copy()
		layer.coeffs = an
		mode_layers.append(layer)
		return mode_layers

DefaultSolver = LayeredSolver
