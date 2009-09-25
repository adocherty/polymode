# _*_ coding=utf-8 _*_
#
#Copyright © 2008 Andrew Docherty

#This program is part of ABCSolver.
#ABCSolver is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

'''
Lightwieght flexible coordinate classes for vector calculations

'''
import logging

from numpy import *
from cached_calculation import *
from ..difflounge import finitedifference

def vector_polar_to_cartesian(At, coord):
    rm,phim = coord.polar2d()
    At = asarray(At)
    A0 = At[0]*cos(phim) - At[1]*sin(phim)
    At[1] = At[0]*sin(phim) + At[1]*cos(phim)
    At[0] = A0
    return At

class Coord(object):
    polar_coordinate = False

    def copy(self, **changed):
        "Create new WgCoord object with parameters of current coord object"
        #Copy self to new coord
        nc = self.__class__()
        
        #Copy old paramters
        nc.__dict__.update(self.__dict__)
            
        #Update changed parameters
        nc.__dict__.update(changed)
        return nc       

    def dot(self, A, B):
        """
        A component dot product of the vectors A and B over the first index
        """
        return (asarray(A)*asarray(B)).sum(axis=0)

    def cross(self, A, B):
        """
        A component cross product of the vectors A and B over the first index
        """
        return cross(A,B,axis=0)

    def fourier_resample(self, data, dcoord, ext=(None,None),
                                rlimits=None, m0=0, fourier=False):
        "Resample fourier data and coord pair to this coord"
        #Get r, phi pairs for all points
        rm, phim = self.polar2d()
        data = asarray(data)
        
        #Sort by rm - we will interpolate the data for each r only once
        sri = argsort(rm.ravel())
        Ncoo = len(rm.ravel())
        if data.ndim>2:
            Ndata = data.shape[0]
            dshape = (Ndata,) + shape(rm)
        else:
            Ndata = 1
            dshape = shape(rm)
        dsamp = zeros((Ndata,Ncoo), complex_)
        mav = zeros(data.shape[-1], complex_)
        
        #Right and left resampling for the coordinate
        if hasattr(dcoord, 'ms'):
            ms = dcoord.ms + m0
            Naz = dcoord.Naz
        else:
            ms = 2*pi*dcoord.wv + m0
            Naz = dcoord.shape[1]

        #Phi shift
        phishift = 0 #pi/dcoord.symmetry

        #Get min/max r for extensions or set and min/max of the calculated field
        if rlimits is not None:
            rmin, rmax = rlimits
        elif hasattr(ext[0], 'rc') and hasattr(ext[1], 'rc'):
            rmin, rmax = ext[0].rc, ext[1].rc
        else:
            rmin, rmax = dcoord.rv[0], dcoord.rv[-1]

        #Use the vector method of the ext instead of their value.
        #Note: evaluating the ext fields node by node here slows things down a lot!
        use_ext_vec = hasattr(ext[0], 'vector') and hasattr(ext[1], 'vector')
        use_ext = (ext[0] is not None) and (ext[1] is not None)

        #Loop over r
        number_different=0; last_r = 0
        for ii in range(Ncoo):
            sii = sri[ii]; r = rm.flat[sii]; phi = phim.flat[sii]

            #DFT coefficients at phi
            fcoeffs = exp(1j*ms*phi)/Naz

            #interpolate using the interior extension if r less than the min
            if r<=rmin:
                if use_ext_vec:
                    mav = ext[0].vector(r)[:Ndata]
                elif use_ext:
                    mav = ext[0]; fcoeffs = 1
                else:
                    mav=asarray(0); fcoeffs = 1

            #interpolate using the exterior extension if r greater than the max
            elif r>=rmax:
                if use_ext_vec:
                    mav = ext[1].vector(r)[:Ndata]
                elif use_ext:
                    mav = ext[1]; fcoeffs = 1
                else:
                    mav=asarray(0); fcoeffs = 1

            #Only sample within calculation domain
            elif abs(r-last_r)>1e-12:
                last_r = r; number_different += 1

                #Find closest r in mode coord
                rii = nonzero(r>dcoord.rv)[0][-1]
                r1,r2 = dcoord.rv[rii:rii+2]
                
                #Linear interpolation
                mav = (data[...,rii+1,:]*(r-r1) + data[...,rii,:]*(r2-r))/(r2-r1)

            dsamp[:,sii] = sum(mav*fcoeffs, axis=-1)

        #Restore the shape to the input shape
        dsamp = dsamp.reshape(dshape)
    
        #Make sure that we have a polar coordinate before taking the fft!
        if fourier:
            dsamp = fft.fft(fft.ifftshift(dsamp, axes=[-1]), axis=-1)
        return dsamp

    @property
    def characteristic_length(self):
        "Characteristic cartesian length"
        return (1, 1)

    def bounds(self):
        "Return maximum and minimum bounds of the coordinate"
        pass


class PolarCoord(Coord):
    polar_coordinate = True

    def __init__(self, rrange=(0,1), arange=(-pi,pi), N=None, border=0, difforder=3, dtype=float_):
        if isscalar(N):
            self.shape = (N,N)
        else:
            self.shape = N
        self.rrange = dtype(rrange)
        self.arange = dtype(arange)
        self.border = border
        self.dtype = dtype
        self.difforder = difforder
    
    def __str__(self):
        return "<%s: r=%s, ϕ=%s, N=%s>" % (type(self).__name__, self.rrange, self.arange, self.shape)
    
    def calc_dr(self):
        dr = (self.rrange[1]-self.rrange[0])/self.shape[0]
        return dr
    dr = CachedCalculation(calc_dr)
    
    def calc_dphi(self):
        dphi = (self.arange[1]-self.arange[0])/self.shape[1]
        return dphi
    dphi = CachedCalculation(calc_dphi)
    
    def calc_rv(self, border=None):
        if border is None: border=self.border
        Nr = self.shape[0]
        rv = arange(border,Nr-border)*self.dr + self.rrange[0]
        return rv
    rv = CachedCalculation(calc_rv)

    def calc_rv_interval(self, border=None):
        if border is None: border=self.border
        Nr = self.shape[0]
        rv = (arange(border,Nr+1-border)-0.5)*self.dr + self.rrange[0]
        if border==0: rv[0]=0
        return rv
    rv_interval = CachedCalculation(calc_rv_interval)

    def calc_phiv(self):
        Naz = self.shape[1]
        phiv = arange(-(Naz//2),(Naz+1)//2)*self.dphi
        return phiv
    phiv = CachedCalculation(calc_phiv)

    def calc_phiv_interval(self):
        Naz = self.shape[1]
        phiv = (arange(-(Naz//2),(Naz+1)//2+1)-0.5)*self.dphi
        return phiv
    phiv_interval = CachedCalculation(calc_phiv_interval)

    def calc_wv(self):
        "Calculate Fourier domain frequencies"
        Naz = self.shape[1]
        wv = arange(-(Naz//2),(Naz+1)//2)/(Naz*self.dphi)
        return fft.fftshift(wv)
    wv = CachedCalculation(calc_wv)

    def bounds(self):
        "Return maximum and minimum bounds of the coordinate"
        return append(self.rrange, self.arange)

    def polar_bounds(self):
        "Return maximum and minimum bounds in terms polar coordinates"
        return append(self.rrange, self.arange)

    # +-----------------------------------------------------------------------+
    # | Conversion/Matrix functions
    # +-----------------------------------------------------------------------+
    def polar2d(self, interval=False):
        if interval:
            rm = self.rv_interval[:,newaxis] + 0*self.phiv_interval[newaxis,:]
            phim = 0*self.rv_interval[:,newaxis] + self.phiv_interval[newaxis,:]
        else:
            rm = self.rv[:,newaxis] + 0*self.phiv[newaxis,:]
            phim = 0*self.rv[:,newaxis] + self.phiv[newaxis,:]
        return rm,phim
        
    def cartesian2d(self, interval=False):
        rm,phim = self.polar2d(interval)
        xm,ym = rm*cos(phim), rm*sin(phim)
        return xm,ym
        
    def convert_transverse_vector(self, At):
        return At
    
    @property
    def characteristic_length(self):
        "Characteristic cartesian length"
        return (self.dr, self.dphi*max(self.rrange))
    
    # +-----------------------------------------------------------------------+
    # | Vector calculus
    # +-----------------------------------------------------------------------+
    def calc_diffr(self):
        D = finitedifference.DifferenceMatrix(2, bandwidth=self.difforder, X=self.rv, dtype=complex)
        D.generate()
        return D
    diffr = CachedCalculation(calc_diffr)

    def calc_diffphi(self):
        D = finitedifference.DifferenceMatrix(2, bandwidth=self.difforder, X=self.phiv, dtype=complex)
        D.generate()
        return D
    diffphi = CachedCalculation(calc_diffphi)
    
    def grad_t(self, f, m0=0, fourier=False):
        rm = self.rv[:,newaxis]
        Gr = self.diffr.diff1(f)
        Gphi = self.diffphi.diff(f, axis=1)/rm
        if self.border==0: Gr[0] = Gr[1]
        return Gr, Gphi
        
    def div_t(self, At, m0=0, fourier=False):
        Ar, Aphi = At
        rm = self.rv[:,newaxis]
        D = (self.diffr.diff1(rm*Ar) + self.diffphi.diff(Aphi, axis=1))/rm
        return D
        
    def curl_t(self, At, m0=0, fourier=False):
        Ar, Aphi = At
        rm = self.rv[:,newaxis]
        Cz = self.diffr.diff1(rm*Aphi)/rm - self.diffphi.diff(Ar, axis=1)/rm
        return Cz

    def int_dA(self, f, m0=0, fourier=False):
        "Integrate in polar coords over an area"
        Ia = (self.rv[:,newaxis]*f).sum()*self.dr*self.dphi
        return Ia

class CartesianCoord(Coord):
    """
    Cartisian coordinates in range
    x = -X ... X
    y = -Y ... Y
    with N points
    """
    polar_coordinate = False

    def __init__(self, X=(-1,1), Y=(-1,1), N=None, difforder=5, border=0, dtype=float_):
        if isscalar(N):
            self.shape = (N,N)
        else:
            self.shape = N
        if isscalar(X):
            self.rangex = dtype([-X,X])
        else:
            self.rangex = dtype(X)
        if isscalar(Y):
            self.rangey = dtype([-Y,Y])
        else:
            self.rangey = dtype(Y)
        
        self.difforder = difforder
        self.border = border
        self.dtype = dtype

    def __str__(self):
        return "<%s: x=%s, y=%s, N=%s>" % (type(self).__name__, self.rangex, self.rangey, self.shape)

    def calc_dx(self):
        dr = (self.rangex[1]-self.rangey[0])/(self.shape[0]-1)
        return dr
    dx = CachedCalculation(calc_dx)

    def calc_dy(self):
        dr = (self.rangey[1]-self.rangey[0])/(self.shape[1]-1)
        return dr
    dy = CachedCalculation(calc_dy)
    
    def calc_xv(self, border=None):
        Nx = self.shape[0]
        if border is None: border=self.border

        xv = arange(border,Nx-border)*self.dx + self.rangex[0]
        return xv
    xv = CachedCalculation(calc_xv)
    
    def calc_yv(self, border=None):
        Ny = self.shape[1]
        if border is None: border=self.border

        xv = arange(border,Ny-border)*self.dy + self.rangey[0]
        return xv
    yv = CachedCalculation(calc_yv)

    def calc_xinterval(self, border=None):
        Nx = self.shape[0]
        if border is None: border=self.border
        xv = (arange(border,Nx+1-border) - 0.5)*self.dx + self.rangex[0]
        return xv
    xinterval = CachedCalculation(calc_xinterval)
    
    def calc_yinterval(self, border=None):
        Ny = self.shape[1]
        if border is None: border=self.border
        xv = (arange(border,Ny+1-border) - 0.5)*self.dy + self.rangey[0]
        return xv
    yinterval = CachedCalculation(calc_yinterval)

    def bounds(self):
        "Return maximum and minimum bounds of the coordinate"
        return append(self.rangex, self.rangey)

    def polar_bounds(self):
        "Return maximum and minimum bounds in terms polar coordinates"
        rs = [hypot(x,y) for x in self.rangex for y in self.rangey]
        rrange = [min(rs), max(rs)]
        
        #If the bounds don't encircle the orgin calculate the max/min args
        if 0:
            args = [arctan2(y,x) for x in self.rangex for y in self.rangey]
            arange = [min(args), max(args)]
        else:
            arange = [-pi,pi]
        return append(rrange, arange)

    # +-----------------------------------------------------------------------+
    # | Conversion/Matrix functions
    # +-----------------------------------------------------------------------+
    def polar2d(self, interval=False):
        xm,ym = self.cartesian2d(interval=interval)
        rm = hypot(xm,ym)
        phim = arctan2(ym,xm)
        return rm,phim
    
    def cartesian2d(self, interval=False):
        if interval:
            xm,ym = self.xinterval[:,newaxis], self.yinterval[newaxis,:]
        else:
            xm,ym = self.xv[:,newaxis], self.yv[newaxis,:]
        return xm+0*ym, 0*xm+ym

    def convert_transverse_vector(self, At):
        rm,phim = self.polar2d()
        return At[0]*cos(phim) - At[1]*sin(phim), At[0]*sin(phim) + At[1]*cos(phim)

    @property
    def characteristic_length(self):
        "Characteristic cartesian length"
        return (self.dx, self.dy)
    
    # +-----------------------------------------------------------------------+
    # | Vector calculus
    # +-----------------------------------------------------------------------+
    def calc_diffx(self):
        D = finitedifference.DifferenceMatrix(2, bandwidth=self.difforder, X=self.xv, dtype=complex)
        D.generate()
        return D
    diffx = CachedCalculation(calc_diffx)

    def calc_diffy(self):
        D = finitedifference.DifferenceMatrix(2, bandwidth=self.difforder, X=self.yv, dtype=complex)
        D.generate()
        return D
    diffy = CachedCalculation(calc_diffy)
        
    def grad_t(self, f, m0=0, fourier=False):
        "Transverse gradiant = ∇⟂f"
        if fourier: raise NotImplementedError
        return (self.diffx.diff1(f), self.diffy.diff(f, axis=1))
    def div_t(self, At, m0=0, fourier=False):
        "Transverse divergence = ∇⟂∙A⟂"
        if fourier: raise NotImplementedError
        return self.diffx.diff1(At[0]) + self.diffy.diff(At[1], axis=1)
    def curl_t(self, At, m0=0, fourier=False):
        "z-component of transverse curl = ẑ∙∇⟂×A⟂"
        if fourier: raise NotImplementedError
        return self.diffx.diff1(At[1]) - self.diffy.diff(At[0], axis=1)

    def int_dA(self, f):
        "Integrate in over an area"
        Ia = f.sum()*self.dx*self.dy
        return Ia

class CartesianCoordZero(CartesianCoord):
    """
    Cartisian coordinates in range
    x = -X ... X-dx
    y = -Y ... Y-dy
    with N points
    """
    def calc_dx(self):
        dr = (self.rangex[1]-self.rangey[0])/self.shape[0]
        return dr
    dx = CachedCalculation(calc_dx)

    def calc_dy(self):
        dr = (self.rangey[1]-self.rangey[0])/self.shape[1]
        return dr
    dy = CachedCalculation(calc_dy)



    
if __name__=="__main__":
    cc = CartesianCoord( rangex=(0,1), rangey=(0,1), N=(5,5), border=0)
    print cc
    print cc.xv
    print cc.yv

    cc2 = CartesianCoord( rangex=(2,4), rangey=(-2,0), N=(5,5), border=0)
    print cc2
    print cc2.xv
    print cc2.yv

