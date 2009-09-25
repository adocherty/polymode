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

"""
Waveguide class to contain all information on the waveguide structure.
Each waveguide had a material, an exterior material and one or more shapes within it.

"""

from __future__ import division
import logging

import numpy as np
from numpy import *

from .difflounge import finitedifference, boundary
from .mathlink import fftp, misc, coordinates, utf8out
from . import Material

class WgCoord(coordinates.Coord):
    '''
    kwargs are:
    - Nshape: number of points in radial & azimuthal direction
    - border:
    - oversample:
    - m0: symmetry class (only used for differentiation)
    - bandwidth: radial finite diff bandwidth (only used for differentiation)
    '''
    polar_coordinate = True
    precalc_names = ['dr','dphi','dA','rv','ms','phiv','rv_interval','phiv_interval','shape','diffr','phirange','rrange','arange']
    rmin = 0
    def __init__(self, rmax, rmin=0, symmetry=1, sectors=1, reflect=0, Nshape=None, border=0, bandwidth=3):
        #self.shape = Nshape
        self.rmax,self.rmin = rmax,rmin
        self.symmetry = symmetry
        self.border = border
        self.bandwidth = bandwidth

        if Nshape is not None: self.Nr,self.Naz = Nshape

        #Modifiers
        self.sectors = sectors
        self.reflect = reflect
        self.phi_extra = 0
        self.shift = 0
        self.os_r = self.os_az = 1

    def set_oversample(self, oversample):
        if size(oversample)==2:
            self.os_r, self.os_az = oversample
        else:
            self.os_r = self.os_az = oversample         

    def clear_oversample(self):
        self.os_r = self.os_az = 1

    def new(self, rmax=None, rmin=None, symmetry=None, Nr=None, Naz=None,
            Nshape=None, sectors=None, reflect=None, border=None):
        "Create new WgCoord object with parameters of current coord object"
        if rmax is None:                rmax = self.rmax
        if rmin is None:                rmin = self.rmin
        if symmetry is None:    symmetry = self.symmetry
        if reflect is None:         reflect = self.reflect
        if sectors is None:         sectors = self.sectors
        if border is None:          border = self.border

        if Nshape is None:          Nshape = (self.Nr,self.Naz)
        if Nr is not None:          Nshape = (Nr, Nshape[1])
        if Naz is not None:         Nshape = (Nshape[0], Naz)

        return WgCoord(rmax,rmin, symmetry, sectors, reflect, Nshape, border)

    # +-----------------------------------------------------------------------+
    # | Pickling marshalling functions
    # | We don't pickle the cached data to give compact storage
    # +-----------------------------------------------------------------------+

    def __getstate__(self):
        "Pickle all needed data, ignore cached data"
        state = self.__dict__.copy()
        #Remove cached information
        for key in state.keys():
            if key.startswith('_'):
                del state[key]
        return state
    
    def __setstate__(self,state):
        "Restore needed data"
        self.__dict__.update(state)
        
    def __hash__(self):
        #Note we're not checking for modifiers like oversampling etc yet ..
        #Currently oversampling (through coord) is depreciated
        return hash((self.rmax, self.symmetry, self.border, self.Nr, self.Naz))
    
    # +-----------------------------------------------------------------------+
    # | Cacheable calculation functions
    # +-----------------------------------------------------------------------+

    def calc_shape(self):
        shape = (self.Nr*self.os_r+2*(1-self.border), self.Naz*self.os_az)
        return shape

    def calc_dr(self):
        dr = (self.rmax-self.rmin)/(self.Nr*self.os_r)
        return dr
    
    def calc_dphi(self):
        if self.reflect:
            dphi = pi/self.symmetry/( ((self.Naz*2-1)*self.os_az)//2)
        else:
            dphi = 2*pi/self.symmetry/(self.Naz*self.os_az)
        return dphi
    
    def calc_dA(self):
        return self.dr*self.dphi
    
    def calc_rv(self, border=None):
        if border is None:
            border=self.border
        rv = arange(border,self.Nr*self.os_r+2-border)*self.dr + self.rmin
        return rv

    def calc_rv_interval(self, border=None):
        if border is None: border=self.border
        rv = (arange(border,self.Nr*self.os_r+1+2-border) - 0.0)*self.dr + self.rmin

        #if border==0: rv[0]=0
        return rv

    def calc_phirange(self):
        Naz = self.Naz
        phirange = -self.dphi*(Naz//2),self.dphi*(Naz//2)
        return phirange

    def calc_rrange(self):
        rrange = (self.border*self.dr+self.rmin,
                    (self.Nr*self.os_r+2-self.border)*self.dr+self.rmin)
        return rrange

    def calc_arange(self):
        arange = [-pi*self.sectors/self.symmetry,pi*self.sectors/self.symmetry]
        return arange

    def calc_ms(self, sectors=None):
        #Calculate the m mu fourier base
        if sectors is None: sectors=self.sectors
        Nazo = sectors*self.shape[1]
        
        #If we have reflection symmetry only calculate positive fouier components
        if self.reflect:
            ms = arange(0,Nazo)*self.symmetry/sectors
        else:
            symmetric = False
            if symmetric and self.Naz%2==0:
                ms = r_[0:Nazo//2,0,-Nazo//2+1:0]*self.symmetry/sectors
            else:
                ms = r_[0:(Nazo+1)//2,-(Nazo//2):0]*self.symmetry/sectors
            
        return ms

    def calc_phiv(self, sectors=None):
        if sectors is None: sectors=self.sectors
        Nazo = self.shape[1]

        #self._phiv = r_[0:(Nazo+1)//2,-(Nazo//2):0]*self.dphi
        if self.reflect:
            phiv = arange(0, Nazo)*self.dphi
        else:
            phiv = arange(-(Nazo//2)-Nazo*((sectors-1)//2), \
                (Nazo+1)//2+Nazo*(sectors//2)+self.phi_extra)*self.dphi
        return phiv + self.shift

    def calc_phiv_interval(self, sectors=None):
        if sectors is None: sectors=self.sectors
        Nazo = self.shape[1]

        phiv = (arange(-(Nazo//2)-Nazo*((sectors-1)//2), \
                (Nazo+1)//2+Nazo*(sectors//2)+1+self.phi_extra) - 0.5)*self.dphi
        return phiv + self.shift

    # +-----------------------------------------------------------------------+
    # | Caching control functions
    # +-----------------------------------------------------------------------+

    #Intercept calls for data & construct them if needed
    def __getattr__(self,name):
        if name in WgCoord.precalc_names:
            try:
                #Try to get the data directly
                return object.__getattribute__(self,"_"+name)
            except AttributeError:
                #If not we need to recalculate it
                self.__setattr__("_"+name,object.__getattribute__(self,"calc_"+name)())
                return object.__getattribute__(self,"_"+name)
        else:
            return object.__getattribute__(self,name)

    #Just ensure we don't set anything we shouldn't
    def __setattr__(self,name,value):
        if name in WgCoord.precalc_names:
            logging.warning( "WgCoord: Precalc items cannot be set directly" )
        else:
            #Need to recalculate precalc items when a parameter has changes
            self.__clear__()
            object.__setattr__(self,name,value)
    
    def __clear__(self):
        "Clean up all precalculated data"
        for name in WgCoord.precalc_names:
            if hasattr(self,"_"+name):
                object.__delattr__(self,"_"+name)

    # +-----------------------------------------------------------------------+
    # | Non-cached functions
    # +-----------------------------------------------------------------------+
    def polar_bounds(self):
        return append(self.rrange, self.arange)

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

    @property
    def characteristic_length(self):
        "Characteristic cartesian length"
        return (self.dr, self.dphi*max(self.rrange))
    # +-----------------------------------------------------------------------+
    # | Functions for finding the vector derivatives
    # +-----------------------------------------------------------------------+

    def calc_diffr(self):
        bw = self.bandwidth
        self._diffr = finitedifference.DifferenceMatrix(2, bandwidth=bw, X=self.rv, dtype=complex_)
        return self._diffr.diff
    
    def diffphi(self,f, n=1, axis=1, m0=0, fourier=False):
        "We need a swapaxis here to correctly multiply m by f if axis<>-1"
        if fourier:
            return (1j*(self.ms+m0))**n * f
        else:
            df = (1j*(self.ms+m0))**n * fftp.fft(f*exp(-1j*m0*self.phiv),axis=axis)
            return fftp.ifft(df*exp(1j*m0*self.phiv),axis=axis)

    def grad_t(self, f, m0=0, fourier=False):
        "Calculate the transverse gradient of f in cylindrical polars"
        rm = self.rv[:,newaxis]

        Gr = self.diffr(f,1)
        Gphi = self.diffphi(f, m0=m0, fourier=fourier)/rm
        if self.border==0: Gr[0] = Gr[1]
        return Gr, Gphi

    def div_t(self, At, m0=0, fourier=False):
        "Calculate the transverse divergence for At in cylindrical polars"
        rm = self.rv[:,newaxis]

        D = (self.diffr(rm*At[0],1) + self.diffphi(At[1], m0=m0, fourier=fourier))/rm
        if self.border==0: D[0] = D[1]
        return D
        
    def curl_t(self, At, m0=0, fourier=False):
        "Calculate the z component of the transverse curl for At in cylindrical polars"
        rm = self.rv[:,newaxis]

        Cz = self.diffr(rm*At[1],1)/rm - self.diffphi(At[0], m0=m0, fourier=fourier)/rm
        if self.border==0: Cz[0] = Cz[1]
        return Cz

    def dot_t(self, A, B):
        "A component dot product of the vectors A and B over the first index"
        return (asarray(A)*asarray(B)).sum(axis=0)
    dot = dot_t
    
    def cross_t(self, A, B):
        "A component cross product of the vectors A and B over the first index"
        return cross(A,B,axis=0)
    cross = cross_t
    
    def int_dA(self, f, full=1):
        "Calculate the area integral of f in cylindrical polars"
        Ia = (self.rv[:,newaxis]*f).sum()*self.dr*self.dphi
        if full:
            Ia *= self.symmetry/self.sectors
        return Ia

# +-----------------------------------------------------------------------+
# | Main WG Shape object
# +-----------------------------------------------------------------------+

class WgShape(object):
    plot_resolution = 0.01          #Resolution for polygon segmentation for plotting
    numerical_difference = False
    zorder = 0      
    def __init__(self, material, center=(0,0), rot=0, zorder=0, xy=False):
        self.material = material
        self.zorder = zorder
        
        #If store center in polar coordinates
        self.center = (0,0)
        self.set_center(center, xy=xy)
        
        #Store rotation in cartesian coordinates
        self.rotation = 0
        self.set_rotation(rot, xy=xy)
        
        if not isinstance(material, Material.Material):
            logging.error( "WgShape: Material not recognised!" )
    
    #ID for shape
    def __hash__(self):
        return hash((self.center[0], self.center[1], self.rotation, self.zorder, self.material))
        
    #Sort/Compare by zorder
    def __cmp__(self,other):
        if hasattr(other, 'zorder'):
            return self.zorder-other.zorder
        else:
            return self.zorder-other

    def __str__(self):
        return "%s: %s" % (self.name, self.material)
    
    def rotate_cartesian(self, xy, rot):
        "Rotate cartesian coordinate tuple of xy by phi"
        xy_ = (cos(rot)*xy[0] - sin(rot)*xy[1], sin(rot)*xy[0] + cos(rot)*xy[1] )
        return xy_
    
    def center_shift(self, dphi=None):
        "Return the polar center (r0,phi0) of the object"
        r0,phi0 = self.center

        #Calculate azimuthal center as a shift + a remainder
        #phi0 = phi_shift*dphi + phi_remainder
        #This, of couse, only works on equi-spaced azimuthal grids
        if dphi is not None:
            phi_shift = (phi0+dphi/2)//dphi
            phi0 = phi0-phi_shift*dphi
            return (r0,phi0,phi_shift)
        else:
            return (r0,phi0)

    #Helper routines to return Polar center & rotation and coordinates
    def convert_to_cartesian(self, p, xy=False):
        "Convert the coordinate tuple to cartesian coordinates"
        if xy:
            return p
        else:
            return p[0]*cos(p[1]), p[0]*sin(p[1])

    def convert_to_polar(self, p, xy=True):
        "Convert the coordinate tuple to polar coordinates"
        if xy:
            r0 = hypot(p[0],p[1])
            phi0 = arctan2(p[1],p[0])
            return (r0,phi0)
        else:
            return p

    @deprecate
    def cartesian_to_polar(self, xy):
        "Convert cartesian coordinate tuple of xy to polar coords"
        r0 = hypot(xy[0],xy[1])
        phi0 = arctan2(xy[1],xy[0])
        return (r0,phi0)

    @deprecate
    def polar_to_cartesian(self, rp):
        "Convert cartesian coordinate tuple of polar coords to cartesian"
        xy = rp[0]*cos(rp[1]), rp[0]*sin(rp[1])
        return xy

    @deprecate
    def polar_rotation(self, phi0):
        "Return the rotation from the polar radial direction"
        rot = self.rotation - phi0
        return rot

    @deprecate
    def cartesian_rotation(self, phi0):
        "Return the rotation from the cartesian x-direction"
        rot = self.rotation
        return rot

    @deprecate
    def shift_rotate_cartesian(self, coord):
        #Shift
        r0, phi0, phi_shift = self.center_shift( dphi=coord.dphi )
        xm = coord.rv[:,newaxis]*cos(coord.phiv-phi0)-r0
        ym = coord.rv[:,newaxis]*sin(coord.phiv-phi0)
        
        #Rotate
        rot = self.get_rotation(phi0, xy=True)
        xmr = cos(rot)*xm + sin(rot)*ym; ymr = - sin(rot)*xm + cos(rot)*ym
        return xmr, ymr, phi_shift

    #Node utility functions
    def nodes_polar_to_cartesian(self,pnodes):
        xynodes = []
        for n in pnodes:
            xynodes += [[n[0]*cos(n[1]), n[0]*sin(n[1])]]
        return xynodes

    def nodes_cartesian_to_polar(self,xynodes):
        pnodes = []
        for n in xynodes:
            pnodes += [[hypot(n[0],n[1]), arctan2(n[1], n[0])]]
        return pnodes

    def set_index_function(self, fn, background=None):
        """
        Store a function fn(d) where d is the distance to the center of the object
        that defines a relative change in the refractive index
        """
        self.index_function = fn
        if background:
            self.background = background    

    def get_center(self, xy=False):
        if xy:
            return self.convert_to_cartesian(self.center)
        else:
            return self.center

    def set_center(self, center, xy=False):
        if xy:
            self.center = self.cartesian_to_polar(center)
        else:
            self.center = center

    def get_rotation(self, c1=None, xy=True):
        c1 = self.center[1] if c1 is None else c1
        if xy:
            return self.rotation
        else:
            return self.rotation - c1

    def set_rotation(self, rot, xy=True):
        if xy:
            self.rotation = rot
        else:
            self.rotation = rot + self.center[1]

    def transform_center(self, shift=(0,0), rotate=0, scale=(0,0), xy=True):
        "Transform center of object, leaving shape unchanged"
        c = self.get_center(xy=xy)
        c = asarray(c)*asarray(scale)
        
        if xy:
            c += array(shift)
            c = self.cartesian_to_polar(c)
            c = (c[0], c[1]+rotate)
        else:
            c = (c[0], c[1]+rotate)
            c = self.convert_to_cartesian(c)
            c += array(shift)
        
        self.set_center(c, xy=not xy)

    def transform(self, shift=(0,0), rotate=0, scale=(0,0), xy=True):
        "Transform object, leaving center unchanged"
        pass

    #External shape interfaces for actual discretization
    def extents(self, xy=False):
        "Return maximum and minimum radial and azimuthal extents of object"
        pass
    
    def get_plot_object(self,rot,**style):
        pass
    
    def get_annulus(self, phi, dphi):
        "Return annuli centered at phi of width dphi"
        pass

    def get_index(self, coord):
        "Return index sampled on a grid"
        pass

    def get_polygons(self, coord):
        "Return shape as polygons"
        pass

    def get_bitmap(self, coord):
        pass

    def dilate(self, distance):
        "Dilate object by distance"
        pass

    def contract(self, distance):
        "Contract object by distance (opposite to dilate)"
        pass

    def copy(self):
        "Return a clone of the shape"
        pass


# +-----------------------------------------------------------------------+
# | Bitmapped shapes
# +-----------------------------------------------------------------------+

class Image(WgShape):
    '''
    Shape constructed from image.
    '''
    name = "Image"
    numerical_difference = True

    def __init__(self, material, image=None, average_over_sectors=None, border=0):
        WgShape.__init__(self,material,(0,0), xy=False)
        self.image = image
        self.border = border
        self.average_over_sectors = average_over_sectors
    
    def extents(self, xy=False):
        "Return radial and azimuthal extents of object"
        R = self.image.rmax
        if xy:
            bounds = (-R,R,-R,R)
        else:
            bounds = (0,R, 0, 2*pi)
        return bounds

    def get_index(self, coord):
        index_data = self.image.regrid(coord, sectors=self.average_over_sectors, rborder=self.border)
        return index_data
    
    #Create polar bitmapped image
    def get_plot_object(self,rot,**style):
        from matplotlib import patches

        radius = self.image.rmax
        r0,phi0 = self.get_center()
        objrot = self.get_rotation(xy=True) + rot

        center = r0*cos(phi0+rot), r0*sin(phi0+rot)
        c = patches.Ellipse(center,2*radius,2*radius,angle=180*(objrot)/pi, **style)
        return [c]

# +-----------------------------------------------------------------------+
# | Basic waveguide shapes
# +-----------------------------------------------------------------------+
class WgAnnulus(object):
    def __init__(self, shape, r=(0,0), phi=(0,0)):
        self.r = (min(r),max(r))
        self.phi = (min(phi),max(phi))
        self.true_phi = (min(phi),max(phi))

        #Copy some information from shape .. these must be named the same for correct merging
        self.material = shape.material
        self.zorder = shape.zorder
        self.center = shape.center
        
        #Optional parameters
        if hasattr(shape, 'index_function'):
            self.index_function = shape.index_function
        if hasattr(shape, 'background'):
            self.background = shape.background

    def get(self):
        return self.r + self.phi

    def get_rstart(self):
        return self.r[0]
    rstart = property(get_rstart)
    
    def get_rend(self):
        return self.r[1]
    rend = property(get_rend)

    def distance(self, rv):
        '''
        Distance from center of original object to the annulus evaluated at rv
        '''
        #Note that we use true_phi, which isn't wrapped to get the correct
        #phi location of the annulus
        r0,phi0 = self.center
        dist = sqrt(r0**2+rv**2 - 2*r0*rv*cos(phi0-average(self.true_phi)))
        return dist

    def wrap_to_sector(self, pmin, pmax, dphi):
        '''
        Wrap the phi coord of the sector to the range between -plim to plim
        This only changes phi, not true_phi
        '''
        p = (self.phi[1]+self.phi[0])/2.
        wid = (self.phi[1]-self.phi[0])/2.
        
        #Wrap midpoint to the phi limits
        pwrap = mod(p-pmin,pmax-pmin+dphi)+pmin
        self.phi = (pwrap-wid,pwrap+wid)

    def __cmp__(self,other):
        "Compare by start radius"
        return cmp(self.rstart, other.rstart)

    def __str__(self):
        return "(%s, z:%d r:%.2g→%.2g, ϕ:%.2g→%.2g)" \
            % ((self.material,self.zorder) + self.r + self.phi)
    def __repr__(self):
        return "(%s, z:%d r:%.2g→%.2g, ϕ:%.2g→%.2g)" \
            % ((self.material,self.zorder) + self.r + self.phi)


class Circle(WgShape):
    """
        Waveguide Circle
        material: the composition of the circle
        center:   (r0,phi0) is the centre of the circle
        radius:   the circle radius
    """
    name = "Circle"
    def __init__(self, material, center=(0,0), radius=1.0, zorder=0, xy=False):
        WgShape.__init__(self,material,center,xy=xy,zorder=zorder)
        self.radius = radius

    #Return the maximum r/phi limits of the object
    #Problem if self.center[0] == 0!
    def extents(self, xy=False):
        "Return radial and azimuthal extents of object"
        Dr = self.radius
        center =self.get_center(xy=xy)
        if xy:
            bounds = (center[0]-Dr, center[0]+Dr,
                    center[1]-Dr, center[1]+Dr)
        else:
            if abs(center[0])<Dr:
                bounds = (0, center[0]+Dr, -pi, pi)
            else:
                Dphi = arcsin(Dr/center[0])
                bounds = (center[0]-Dr, center[0]+Dr,
                    center[1]-Dphi, center[1]+Dphi)
        
        return bounds

    def copy(self):
        newcircle = self.__class__(self.material, self.center)
        newcircle.radius = self.radius
        newcircle.zorder = self.zorder
        return newcircle

    def dilate(self, distance):
        "Dilate object by distance"
        self.radius += distance
        
    #Create the circle out of fourier annuli
    def get_annuli(self, phiv, dphi):
        "Decimate object into radial annuli"
        ra,rb,phia,phib = self.extents()
        ax = self.radius
        r0,phi0 = self.center
        
        annuli = []
        phiv_inrange = filter(lambda p: phia<=p+dphi/2<=phib, phiv)
        
        if r0<ax:
            for phi in phiv_inrange:
                rb = r0*cos(phi+dphi/2-phi0) + sqrt(ax**2 - r0**2*sin(phi+dphi/2-phi0)**2)
                annuli += [ WgAnnulus(self, (0,rb), (phi,phi+dphi)) ]
        else:
            for phi in phiv_inrange:
                ra = r0*cos(phi+dphi/2-phi0) - sqrt(ax**2 - r0**2*sin(phi+dphi/2-phi0)**2)
                rb = r0*cos(phi+dphi/2-phi0) + sqrt(ax**2 - r0**2*sin(phi+dphi/2-phi0)**2)
                annuli += [ WgAnnulus(self, (ra,rb), (phi,phi+dphi)) ]
        return annuli

    def get_plot_object(self, rot, **style):
        "Create matplotlib object and return, for plotting"
        from matplotlib import patches
        
        r0,phi0 = self.get_center()
        ax = self.radius
        objrot = self.get_rotation(xy=True) + rot

        center = r0*cos(phi0+rot), r0*sin(phi0+rot)
        e = patches.Ellipse(center,2*ax,2*ax,angle=180*(objrot)/pi, **style)
        return [e]

class Polygon(WgShape):
    """
    Annular Polygon
    nodes: list of nodes as [ (r1,phi1), (r2, phi2) ... ]
    xy: in polar coords if False, cartesian if True
    """
    name = "Polygon"
    def __init__(self, material, center=(0,0), nodes=[], zorder=0, xy=False):
        #Convert cartesian nodes to polar
        self.nodes = [self.convert_to_polar(node, xy=xy) for node in nodes]

        WgShape.__init__(self,material,center,zorder=zorder,xy=xy)

    def set_center(self, center, xy=False):
        #The take the cartesian difference of centers
        dxy = asarray(self.convert_to_cartesian(center, xy=xy))
        dxy -= asarray(self.get_center(xy=True))
        self.transform(shift=dxy)

        WgShape.set_center(self,center,xy)

    def set_rotation(self, rot, xy=True):
        #Rotate nodes by difference
        drot = self.get_rotation(xy=xy) - rot
        self.transform(rotate=drot)

        WgShape.set_rotation(self,rot,xy)
        
    def transform(self, shift=(0,0), rotate=0, scale=(1,1), centerxy=None, xy=True):
        "Transform object, leaving center unchanged"
        if centerxy is None:
            centerxy = asarray(self.get_center(xy=True))

        for kk in range(len(self.nodes)):
            #Subtract center to do local transformations
            nodexy = self.convert_to_cartesian(self.nodes[kk]) - centerxy

            if xy:
                #Scale and shift in cartesians
                nodexy = asarray(nodexy)*asarray(scale)
                nodexy += array(shift)

                #Add rot to angle and back to cartesian
                node = self.convert_to_polar(nodexy)
                nodexy = self.convert_to_cartesian([node[0], node[1]+rotate])
            else:
                #Scale in local polars
                node = self.cartesian_to_polar(nodexy)
                node = asarray(node)*asarray(scale)

                #Rotate in local polars
                nodexy = self.convert_to_cartesian([node[0], node[1]+rotate])
                
                #Shift in cartesians
                nodexy += array(shift)
            self.nodes[kk] = self.convert_to_polar(nodexy+centerxy)

    #Return the maximum r/phi limits of the object
    def extents(self, xy=False):
        "Return radial and azimuthal extents of object"
        if xy:
            na = array(self.nodes_polar_to_cartesian(self.nodes))
        else:
            na = array(self.nodes)
        bounds = (na[:,0].min(), na[:,0].max(), na[:,1].min(), na[:,1].max())
        return bounds

    def copy(self):
        newshape = self.__class__(self.material, self.center)
        newshape.nodes = list(self.nodes)
        newshape.zorder = self.zorder
        return newshape

    def dilate(self, d):
        "Dilate object by a distance"
        def node_angle(na, nb):
            return arctan2( na[1]-nb[1], na[0]-nb[0] )
        
        #The direction of the interior of the polygon
        #Currently assumes that the nodes are ordered clockwise
        edirection = 1
        
        #The routine will fail if there are repeated nodes
        if absolute(array(self.nodes[-1])-array(self.nodes[0])).max()<1e-10:
            self.nodes = self.nodes[:-1]            
        
        Nnodes = len(self.nodes)
        xynodes = self.nodes_polar_to_cartesian(self.nodes)
        newnodes = []
        for kk in range(Nnodes):
            previous = mod(kk-1,Nnodes)
            next = mod(kk+1,Nnodes)
    
            #Angles from/to previous/next nodes
            p1=node_angle(xynodes[previous], xynodes[kk])
            p2=node_angle(xynodes[kk], xynodes[next])

            #Calculate dilation node location (in cartesians)
            A = array([[cos(p1),-cos(p2)],[sin(p1),-sin(p2)]])
            Y = array([d*cos(p2+edirection*pi/2)-d*cos(p1+edirection*pi/2),
                d*sin(p2+edirection*pi/2)-d*sin(p1+edirection*pi/2)])
            
            #Dilation nodes
            dns = np.linalg.solve(A,Y)
            xy = [dns[0]*cos(p1)+d*cos(p1+edirection*pi/2),
                dns[0]*sin(p1)+d*sin(p1+edirection*pi/2)]
            
            #Add new node to list
            newnodes += [ (xynodes[kk][0] + xy[0],xynodes[kk][1] + xy[1]) ]
        self.nodes = self.nodes_cartesian_to_polar(newnodes)

    def get_annuli(self, phiv, dphi):
        "Decimate object into radial annuli"
        nodes = self.nodes
        Nn = len(nodes)

        #Calculate the angle < 2pi
        def subtended(x1, x2):
            phid = x1-x2
            if abs(phid)>pi:
                phid = mod(x1,2*pi) - mod(x2,2*pi)
            return phid
        
        def within(x1, x, x2):
            if abs(x1-x2)<=pi and min(x1,x2)<=x<=max(x2,x1):
                return True
            x1=mod(x1,2*pi); x2=mod(x2,2*pi); x=mod(x,2*pi)
            if abs(x1-x2)<=pi and min(x1,x2)<=x<=max(x2,x1):
                return True
            return False
        
        #Check if we enclose the origin
        phispan = 0
        for n in range(Nn):
            r1, phi1 = nodes[n]; r2, phi2 = nodes[(n+1)%Nn]
            phispan += subtended(phi2,phi1)
        
        #If phispan==0, not enclosing, if phispan==2pi, enclosing
        interior=False
        if phispan>pi:
            interior=True
            phia=-pi; phib=pi
            
        #Construct annuli from polygon
        annuli = []
        for phi in phiv:
            #find r intersects:
            risect = []
            for k in range(Nn):
                #is phi in range?
                r1, phi1 = nodes[k]
                r2, phi2 = nodes[(k+1)%Nn]
                
                if within(phi1,phi+dphi/2,phi2):
                    risect += [ r1*r2*sin(subtended(phi1,phi2)) \
                        / (r2*sin(subtended(phi+dphi/2,phi2)) \
                        - r1*sin(subtended(phi+dphi/2,phi1))) ]

            if interior and len(risect)>0:
                annuli+=[WgAnnulus(self, (0,max(risect)), (phi,phi+dphi))]
            elif len(risect)>1:
                annuli+=[WgAnnulus(self, (min(risect),max(risect)), (phi,phi+dphi))]

        return annuli

    def to_nodelist(self, rotate=0, close=False):
        "Get nodelisting for export"
        nodes = [(nr,na+rotate) for nr,na in self.nodes]

        if close:
            nodes += nodes[:1]

        xynodes = self.nodes_polar_to_cartesian(nodes)
        return xynodes

    def get_plot_object(self, rot, **style):
        from matplotlib import patches

        #Get cartesian nodes
        xy = self.to_nodelist(rotate=rot, close=True)
        e = patches.Polygon(array(xy), **style)
        return [e]

class RegularPolygon(Polygon):
    """
    Creates a regular polygon with N sides
    Nsides: Number of sides
    length: Side length
    center: center of shape
    """
    name = "RegularPolygon"
    def __init__(self, material, center=(0,0), rot=0, Nsides=4, length=1.0, zorder=0, xy=False):
        #Enclosing circle radius
        a = 2*pi/Nsides
        r = 0.5*length/sin(a/2)
        
        #Create nodes for polygon without rotation at center 0
        phi = arange(0,Nsides)*a + 0.5*a
        xynodes = zip(r*cos(phi), r*sin(phi))
        self.nodes = self.nodes_cartesian_to_polar(xynodes)
        
        #Final setup: apply rotation and center shift
        WgShape.__init__(self, material, center, rot, xy=xy, zorder=zorder)

class Rectangle(Polygon):
    """
    Rectangle(material, center=(0,0), axes=(1.0,1.0), rot=0, xy=False)
    xy: in polar coords if False, cartesian if True
    """
    name = "Rectangle"
    def __init__(self, material, center=(0,0), rot=0, axes=(1.0,1.0), zorder=0, xy=False):
        #Setup some parameters
        self.axes = axes

        #Create nodes for rectangle, without rotation
        ax,ay = axes
        xynodes = [ (-ax/2,-ay/2), (ax/2,-ay/2), (ax/2,ay/2), (-ax/2,ay/2) ]
        self.nodes = self.nodes_cartesian_to_polar(xynodes)
        
        #Final setup: apply rotation and center shift
        WgShape.__init__(self, material, center, rot, xy=xy, zorder=zorder)

class Curve(Polygon):
    """
    Nurbs based curve shape

    Curve(material, nodes, center=(0,0), rot=0, sampling=100, xy=False)
    nodes: list of nodal coordinates as (x,y) pairs
      Must describe a closed shape, although this is not checked.
    center: (r0,phi0) or (x0,y0) coordinate pair. The curve origin is moved to this point.
    xy: use cartesian coords for center if True, polar coords if False
    """
    name = "Nurbs Curve"
    def __init__(self, material, nodes, center=(0,0), rot=0, sampling=100, zorder=0, xy=False):
        #Nodes for shape and knots
        if not xy:
            nodes = self.nodes_polar_to_cartesian(nodes)

        nodes = asarray(nodes).T
        knots = linspace(0.,1.,nodes.shape[1]-1)
        knots = r_[0,0,knots,1,1]

        #Create NURBS curve
        from Nurbs import Crv
        self.curve = Crv.Crv(nodes,knots)

        #Create nodes by sampling NURBS curve, ie convert curve to line segments
        lsnodes=[]
        xynodes = self.curve.pnt3D(linspace(0., 1, sampling)).T
        self.nodes = self.nodes_cartesian_to_polar(xynodes)
        
        #Final setup: apply rotation and center shift
        WgShape.__init__(self, material, center, rot, xy=xy, zorder=zorder)

class Ellipse(Polygon):
    """
        Nurbs based ellipse shape
        axes: (a_r,a_phi) are the ellipse semi-major & semi-minor axes
            (respective to the tangents to phi & r at the ellipse centre)
        rot: is the rotation in radians from the radial direction
        xy: use cartesian coords if True, polar coords if False
    """
    name = "Nurbs Curve"
    def __init__(self, material, center=(0,0), rot=0, axes=(1,1), sampling=50, zorder=0, xy=False):
        #Create NURBS circle
        from Nurbs import Crv
        from Nurbs.Util import translate, scale, rotz
        self.curve = Crv.UnitCircle()

        #Elipse from NURBS circle by a transform
        self.curve.trans(scale(axes))
        
        #Create nodes by sampling NURBS curve, ie convert curve to line segments
        xyznodes = self.curve.pnt3D(linspace(0.,1.,sampling))
        xynodes = transpose(xyznodes[:2])
        
        self.nodes = self.nodes_cartesian_to_polar(xynodes)
        
        #Final setup: apply rotation and center shift
        WgShape.__init__(self, material, center, rot, xy=xy, zorder=zorder)

class ChannelAnnulus(Polygon):
    """
        Annulus with channel of width d at either extent
        r: (r_start,r_end) are the radial extents
        phi: (phi_start, phi_end) are the azimuthal extents
        d1,d1: cartesian distance between phi_start, phi_end & actual annulus
        ni: the refractive index
    """
    name = "Annulus with Channel"
    def __init__(self, material, r=(0,1.0), phi=(1,0), sampling=50, d=0, zorder=0):
        center = (average(r), average(phi))
        WgShape.__init__(self, material, center, zorder=zorder)

        if d==0 or r[0]==0 or r[1]==0:
            dr1a = dr2a = dr1b = dr2b = 0
        elif iterable(d):
            dr1a = arcsin(d[0]/r[0]); dr2a = arcsin(d[0]/r[1])
            dr1b = arcsin(d[1]/r[0]); dr2b = arcsin(d[1]/r[1])
        else:
            dr1a = dr1b = arcsin(d/2/r[0]); dr2a = dr2b = arcsin(d/2/r[1])
        
        self.nodes = []
        dphi = (phi[1]-phi[0]-dr1a-dr1b)/sampling
        for pi in arange(0,sampling+1):
            self.nodes += [ (r[0],phi[0]+dr1a+pi*dphi) ]
        
        dphi = (phi[1]-phi[0]-dr2a-dr2b)/sampling
        for pi in arange(0,sampling+1):
            self.nodes += [ (r[1],phi[1]-dr2b-pi*dphi) ]

class Annulus(WgShape):
    """
        Annulus
        r: (r_start,r_end) are the radial extents
        phi: (phi_start, phi_end) are the azimuthal extents
        ni: the refractive index
    """
    name = "Annulus"
    def __init__(self, material, r=(0,1), phi=(-pi,pi), zorder=0):
        center = (average(r), average(phi))
        WgShape.__init__(self,material,center,xy=False, zorder=zorder)
        self.r = list(r)
        self.phi = list(phi)

    #Return the maximum r/phi limits of the object
    #Problem if self.center[0] == 0!
    def extents(self, xy=False):
        "Return radial and azimuthal extents of object"
        if xy:
            raise NotImplementedError
        else:
            bounds = self.r + self.phi
        return bounds

    def dilate(self, distance):
        "Dilate object by distance"
        self.r[0] -= distance
        self.r[1] += distance
        
    #Create the circle out of fourier annuli
    def get_annuli(self, phiv, dphi):
        "Decimate object into radial annuli"
        ra,rb,phia,phib = self.extents()
        
        annuli = []
        for phi in filter(lambda p: phia<=p+dphi/2<=phib, phiv):
            annuli += [ WgAnnulus(self, (ra,rb), (phi,phi+dphi)) ]
        return annuli

    def get_plot_object(self, rot, **style):
        from matplotlib import patches
        ra,rb,phia,phib = self.extents()
        
        #Draw with min resolution of 0.01 pi
        Nseg = max(ceil((phib-phia)/(0.01*pi)), 10)

        phiv = arange(0,Nseg)*(phib-phia)/(Nseg-1) + phia
        phi = append(phiv, phiv[::-1])
        r = [ra]*Nseg + [rb]*Nseg
        
        #Transform object coordinates
        xy = [r*cos(phi+rot), r*sin(phi+rot)]
        e = patches.Polygon(array(xy).T, **style)
        return [e]

# ---------------------------------------------------------------------
#
#  Shape utility functions
#
# ---------------------------------------------------------------------

def copy_and_rotate_shapes(shapes, rot):
    "Create new shapes with centers rotated by rot"
    newshapes = []
    for shape in atleast_1d(shapes).flat:
        ns = shape.copy()
        c = ns.get_center()
        ns.set_center((c[0], c[1]+rot))
        newshapes.append(ns)
    return newshapes

def transform_shape_centers(shapes, **kwargs):
    '''
    Transform all shapes in list by the following transforms:
    rotate: radians to rotate about origin
    scale: scale center x,y or r,phi by this tuple individually
    shift: shift center x,y or r,phi by this tuple individually
    xy: interpret shift and scale as cartesian if True
    
    Note: None of these transforms change the shape itself,
        only it's center.
    '''
    for shape in atleast_1d(shapes).flat:
        shape.transform_center(**kwargs)
    return shapes

def create_coated_shape(shape, coating_material, coating_thickness, outer=True, inner=False):
    """
    Return a coating to the shape from a specific material
    with specified thickness.

    outer: put coating solely on outside of the shape
    inner: put coating solely on inside of the shape
    
    Note if both outer and inner are true, coating is centered on the shape boundary
    """
    coating = shape.copy()
    coating.material = coating_material
    coating.zorder = shape.zorder-1

    #Contract shape and dilate coating by the appropriate amounts
    if outer and inner:
        coating.dilate(coating_thickness/2)
        shape.dilate(-coating_thickness/2)
    elif outer:
        coating.dilate(coating_thickness)
    elif inner:
        shape.dilate(-coating_thickness)

    return [shape, coating]

def create_square_tiling(shape, offset=1, layers=2, Dx=1, Dy=1, symmetry=4):
    '''
    Create a square lattice with specified layers and spacing
    layers:         number of layers
    offset:         number of layers to skip
                                (default 1 is to remove the central hole)
    Dx:                     x spacing
    Dy:                     y spacing
    symmetry: the waveguide symmetry (1, 2 or 4)
    '''
    assert symmetry in [1,2,4], "Symmetry can only be 1,2 or 4"

    shapes = []

    #All other shapes
    for nr in range(1,layers+offset):
        istart = offset if nr<offset else 0
        for ii in range(istart, layers+offset):
            center = nr*Dx,ii*Dy

            newshape = shape.copy()
            newshape.set_center(center, xy=1)
            shapes.append(newshape)

            if symmetry<4:
                center = -ii*Dx,nr*Dy
                newshape = shape.copy()
                newshape.set_center(center, xy=1)
                shapes.append(newshape)

    if symmetry<2:
        shapes += copy_and_rotate_shapes(shapes, pi)

    #Special case of no central gap
    if offset==0:
        newshape = shape.copy()
        newshape.set_center((0,0))
        shapes.append(newshape)

    return shapes

def create_hexagonal_tiling(shape, offset=1, layers=1, D=1, symmetry=6):
    '''
    Create a hexagonal lattice with specified layers and spacing D
    layers:         number of layers
    offset:         number of layers to skip
                                (default 1 is to remove the central hole)
    Dx:                     x spacing
    Dy:                     y spacing
    symmetry:   waveguide symmetry one of (1,2,3,6)
    '''
    hshapes = []
    #Recursively add tiled shapes from a list of shapes
    if iterable(shape):
        for s in shape:
            hshapes+=create_hexagonal_tiling(s, offset, layers, D, symmetry)
        return hshapes

    #Create hexagonal tiling
    for nl in range(offset,layers+offset):
        for ii in range(nl):
            center = (nl-0.5*ii)*D, ii*sin(pi/3)*D

            newshape = shape.copy()
            newshape.set_center(center, xy=1)
            hshapes.append(newshape)

    #Copy shapes for other symmetries
    if symmetry<6:
        hshapes += copy_and_rotate_shapes(hshapes, pi/3)
    if symmetry<3:
        hshapes += copy_and_rotate_shapes(hshapes, pi/3)
    if symmetry<2:
        hshapes += copy_and_rotate_shapes(hshapes, pi)

    #Special case of no central gap
    if offset==0:
        newshape = shape.copy()
        newshape.set_center((0,0))
        hshapes.append(newshape)

    return hshapes
        
# ---------------------------------------------------------------------
#
#  Classes to combine the shapes in the waveguide in different ways
#
# ---------------------------------------------------------------------
class CombinedIndex(object):
    def __init__(self, coord, mat0, wl, dtype=complex_):
        #Setup the storage area for index caclulations
        nshape = coord.shape
        self.ni2 = empty(nshape, dtype=dtype)
        self.drlogn = zeros(nshape, dtype=dtype)
        self.dazlogn = zeros(nshape, dtype=dtype)

        self.wl = wl
        self.n0 = mat0.index(wl)
        self.ni2[:] = self.n0**2

    def material_index(self, mat):
        return mat.index(self.wl)

    def store_index2(self, x):
        self.ni2[:] += x

    def store_drlogn(self, x):
        self.drlogn[:] += x

    def store_dazlogn(self, x):
        self.dazlogn[:] += x

class CombinedMask(object):
    def __init__(self, coord, mat0, mask_mat, dtype=complex_):
        #Setup the storage area for index caclulations
        nshape = coord.shape
        self.ni2 = zeros(nshape, dtype=dtype)
        self.drlogn = zeros(nshape, dtype=dtype)
        self.dazlogn = zeros(nshape, dtype=dtype)
    
        self.material = mask_mat
        self.n0 = (mat0==self.material)
        self.ni2[:] = self.n0

    def material_index(self, mat):
        return mat==self.material

    def store_index2(self, x):
        self.ni2[:] += x

    def store_drlogn(self, x):
        self.drlogn[:] += x

    def store_dazlogn(self, x):
        self.dazlogn[:] += x

class Combine(object):
    """
    Combine objects take a list of waveguide shapes and combine them
    taking into account the zorder (which shape is behind what) and
    finally generating the index mesh.
    """
    def __init__(self, coord):
        self.coord = coord
    
    def merge_shapes(self, shapes):
        """
        Takes a list of shapes and merges them into the internal list
        """
        pass
    def generate(self, ri):
        """
        Combines the internal shape list with the CombinedIndex object ri
        Stores the index mesh in ri
        """
        pass

class DirectCombine(Combine):
    """
    Combine objects directly returning an index difference profile
    """
    def __init__(self, coord):
        self.coord = coord
        self.shapes = []
        
    def merge_shapes(self, shapes):
        """
        Takes a list of shapes and merges them into the internal list
        """
        self.shapes += shapes
        
    def generate(self, ri):
        """
        Combines the internal shape list with the CombinedIndex object ri
        Stores the index mesh in ri
        """
        for shape in self.shapes:
            ni = ri.material_index(shape.material)
        
            if shape.numerical_difference:
                ni2 = shape.get_index(self.coord)**2*(ni**2 - ri.n0**2)
                ri.store_index2( ni2 )
                ri.store_drlogn( self.coord.diffr(log(ni2+ri.n0**2), 1, axis=0) )
                ri.store_dazlogn( self.coord.diffphi(log(ni2+ri.n0**2), 1, axis=1) )
            else:
                drlogn, dazlogn, ni2 = shape.get_index(self.coord, wl)
                ri.store_index2(ni2*(ni**2 - ri.n0**2))
                ri.store_drlogn(drlogn)
                ri.store_dazlogn(dazlogn)


class AnnularCombine(Combine):
    def __init__(self, coord, symmetry, limit_to_sector=False):
        self.limit_to_sector = limit_to_sector
        self.wrap_to_symmetry = not limit_to_sector
        
        self.symmetry = symmetry
        self.coord = coord
        self.annuli = []

    def annular_overlap(self, a1, a2):
        "Returns true if annuli a1 & a2 overlap"
        if a1.rstart<=a2.rstart<a1.rend:
            return True
        elif a2.rstart<=a1.rstart<a2.rend:
            return True
        return False
        
    def annular_split(self, a1, a2):
        "Splits overlapping annuli into 3 separate annuli"
        #Determine zone 3
        if a1.rend>a2.rend:
            zone3 = a1
        else:
            zone3 = a2

        #Determine zone 2
        if a1.zorder>a2.zorder:
            zone2 = a1
        else:
            zone2 = a2
    
        rendsort = sorted([a1.rend,a2.rend])
    
        #Split into 3 annuli, if they are wider than 0
        na = []
        if a2.rstart>a1.rstart:
            na += [WgAnnulus(a1, (a1.rstart,a2.rstart), a1.phi)]

        na += [WgAnnulus(zone2, (a2.rstart,rendsort[0]), a1.phi)]
    
        if a2.rend<>a1.rend:
            na += [WgAnnulus(zone3, rendsort, a1.phi)]
        return na

    def split_at_phi(self, a1, phi=0):
        "Splits one annulus into two at phi"
        #Split into 3 annuli, if they are wider than 0
        na = [a1]
        if a1.phi[0]<phi<a1.phi[1]:
            na = [WgAnnulus(a1, (a1.rstart,a2.rstart), (a1.phi[0], phi)), \
                   WgAnnulus(a1, (a1.rstart,a2.rstart), (phi, a1.phi[1]))]
        return na

    def annular_merge(self, a1, a2):
        "Merges touching annuli into one"
        na = None
        if a1.rend==a2.rstart and a1.material==a2.material and a1.zorder==a2.zorder:
            na = [WgAnnulus(a1, (a1.rstart,a2.rend), a1.phi)]
        return na
        
    def merge_shapes(self, shapes, oversample=1):
        """
        Takes a list of shapes converts them to annuli and merges these into the internal list
        """
        self.coord.set_oversample(oversample)
        
        #Check limits of phi on objects:
        minphi = min([min(shape.extents()[2:]) for shape in shapes])
        maxphi = max([max(shape.extents()[2:]) for shape in shapes])
    
        #Calculate azimuthal coordinates, calculate annuli over 2*pi range to
        #automatically wrap the shapes to the symmetry sector or just in one
        #sector to discard bits of shapes outside the sector
        dphi = self.coord.dphi
        if self.limit_to_sector:
            phiv = self.coord.calc_phiv()
        else:
            nu = 2*pi/self.coord.symmetry
            maxsector = 2*ceil(max(maxphi/nu,-minphi/nu)-0.5)+1
            phiv = self.coord.calc_phiv(sectors=int(maxsector))
        
        #phirange = self.coord.phirange
        phirange = min(self.coord.phiv), max(self.coord.phiv)
        
        #Collect annuli from all shapes in waveguide 
        new_annuli = self.annuli
        for shape in shapes:
            shape_annuli = shape.get_annuli(phiv, dphi)
    
            ai=0
            while ai<len(shape_annuli):
                a = shape_annuli[ai]

                #Reduce phi in annuli to the symmetry sector
                #This enables out-of-sector objects to wrap and merge correctly
                if self.wrap_to_symmetry:
                    a.wrap_to_sector(phirange[0], phirange[1], dphi)
                
                #If reflection symmetry is on, then be sure to split at zero            
                if self.coord.reflect:
                    shape_annuli[ai:ai+1] = self.split_at_phi(a, 0)

                ai+=1
            
            new_annuli += shape_annuli
        
        #Iterate over all phi and split and merge overlapping annuli
        self.annuli=[]
        for phi in self.coord.calc_phiv():
            #Select annuli at specific phi (assume all the same width)
            pan = filter(lambda a: abs(a.phi[0]-phi)<dphi/10, new_annuli)

            #If more than one annulus at specific phi, split and merge any overlapping
            if len(pan)>1:
                ii=0; pan = sorted(pan)
                while ii<len(pan)-1:
                    #If overlapping split and insert new annuli
                    if self.annular_overlap(pan[ii], pan[ii+1]):
                        nas = self.annular_split(pan[ii], pan[ii+1])
                        pan[ii:ii+2] = nas
                        
                        #If the list has changed so sort and restart
                        #Note: the minimum annulus will remain at ii as
                        # only ii+ are new and all those less than ii are sorted
                        # so we can safely restart at ii, however, decrement by one
                        # so that merging will work
                        ii = max(ii-1,0)
                        pan = sorted(pan)

                    #Otherwise try and merge annuli
                    else:
                        nas = self.annular_merge(pan[ii], pan[ii+1])
                        if nas is not None:
                            pan[ii:ii+2] = nas
                        else:
                            #Here no annuli to split or merge, so move on to the next!
                            ii+=1
            self.annuli += pan

        #Have to remember to clear oversampling!
        self.coord.clear_oversample()
        return self.annuli

    def plot(self, polar=True, color=None):
        import pylab as p_
        from matplotlib import patches
        ax = p_.axes()
        
        rot=0; xmin=ymin=inf; xmax=ymax=0
        for annulus in self.annuli:
            style = {'fill':True, 'fc':annulus.material.color, 'alpha':0.6}
            
            ra, rb, phia, phib = annulus.get()

            #Create the coord pairs for the annulus
            if polar:
                xy=[ (ra*cos(phia), ra*sin(phia)), (ra*cos(phib), ra*sin(phib)), \
                    (rb*cos(phib), rb*sin(phib)), (rb*cos(phia), rb*sin(phia)), \
                    (ra*cos(phia), ra*sin(phia)) ]
            else:
                xy=[ (ra,phia), (ra,phib), (rb,phib), (rb,phia), (ra,phia) ]
            
            xmin = min(xmin, min(array(xy)[:,0]))
            xmax = max(xmax, max(array(xy)[:,0]))
            ymin = min(ymin, min(array(xy)[:,1]))
            ymax = max(ymax, max(array(xy)[:,1]))
            
            #Add a polygon representing the annulus to the list
            ax.add_artist( patches.Polygon(xy, **style) )
            
        ax.axis([xmin, xmax, ymin, ymax])

    def generate(self, ri):
        """
        After all shapes have given their annuli
        * calculate the index with mesh averaging
        * calcualte the differences with mesh averaging
        """
        Nr,Naz = self.coord.rv.shape[0], self.coord.phiv.shape[0]
        dm = self.coord.symmetry; dr = self.coord.dr
        rv = self.coord.rv; ms = self.coord.ms
        n0 = ri.n0
        
        gridav = ones(self.coord.rv.shape)      #Grid averaging weights
        Fth = zeros(self.coord.ms.shape, dtype=complex_)
        for annulus in self.annuli:
            ra, rb, phia, phib = annulus.get()

            #Limit maximum radial extent to be within the maximum radius
            if ra>self.coord.rmax: continue
            if rb>self.coord.rmax: rb = self.coord.rmax
            
            #Get the refractive index of the annulus (via the ri object
            #for flexibility to do other things)
            ni = ri.material_index(annulus.material)
            
            #The radial extent of the annulus in discretized coords
            ka,kb = searchsorted( rv, [ra-dr/2, rb+dr/2] )
            ka = min(ka, Nr-1); kb = min(kb, Nr-1)
            
            #Grid averaging
            gridav[:] = 1
            if ka==kb:              #Annulus is within one grid spacing
                            gridav[ka] = (rb-ra)/dr
            else:                       #Otherwise do both start & ending grid location
                            if ka>0:
                                gridav[ka] = (rv[ka]-(ra-dr/2))/dr
                            if kb<Nr:
                                gridav[kb-1] = ((rb+dr/2)-rv[kb-1])/dr
            
            #Full annulus
            if (phib-phia)<0 or (phib-phia)>=(2*pi/dm):
                iFth = 1; dFth = 0
                
            #Truncated analytic fourier series of step annulus
            #
            # G_m = 1/T ∫_ϕ₀^ϕ₁ g(x) exp(-2πi m/T ϕ) dϕ          with T= 2π/μ
            #     = -μ/(2πi mμ) [exp(-i mμ ϕ₁) - exp(-i mμ ϕ₀)]
            # G_0 = μ/(2π) [ϕ₁ - ϕ₀]
            #
            elif self.coord.reflect:
                raise NotImplementedError, "Reflection symmetry not yet implemented"
            else:
                Fth[1:]=dm/(ms[1:]*2j*pi)*(exp(-1j*ms[1:]*phia)-exp(-1j*ms[1:]*phib))
                Fth[0]=dm*(phib-phia)/(2*pi)

                #Correct for symmetry 
                if self.symmetry!=self.coord.symmetry:
                    Fth *= self.symmetry/self.coord.symmetry
                    set_to_zero = nonzero(mod(ms,self.symmetry)!=0)
                    Fth[set_to_zero] = 0

                #Return aligned to coord phiv
                iFth = fftp.ifft(Fth*Naz)
                dFth = fftp.ifft(1j*ms*Fth*Naz)
                #iFth = fftp.fftshift(fftp.ifft(Fth*Naz))
                #dFth = fftp.fftshift(fftp.ifft(1j*ms*Fth*Naz))
            
            #Add annuli with index function
            if hasattr(annulus,'index_function'):
                #The index of the background
                nb = n0
                if hasattr(annulus, 'background'):
                    nb = ri.material_index(annulus.background)

                #Caclulate refractive index as function of cartesian distance from
                #center of object, set to zero outside object
                #The profile is references to nb not n0 (could be the same)
                dv = annulus.distance(rv[:,newaxis])
                nf2 = gridav[:,newaxis]*(((ni-nb)*annulus.index_function(dv)+nb)**2-n0**2)
                nf2[:ka]=0; nf2[kb:]=0
                nf2 = nf2.astype(complex_)
                
                ri.ni2[ka:kb] += nf2[ka:kb]*iFth
                if hasattr(self,'index_deriv'):
                    ri.drlogn[:] += (2*annulus.index_deriv(dv)/(nf2+n0**2)/(ni**2-n0**2))*iFth
                else:
                    ri.drlogn[ka-1:kb+1] += self.coord.diffr(log((nf2+n0**2)))[ka-1:kb+1]*iFth
                ri.dazlogn[ka:kb] += (log(nf2[ka:kb]+n0**2)-log(n0**2))*dFth

            #Add step annuli with grid averaging
            else:
                ri.ni2[ka:kb] += gridav[ka:kb,newaxis]*(ni**2-n0**2)*iFth
                
                #There are problems with calculating derivatives at the computational
                #boundary so prohibit this for the time being
                if 0<ka<Nr-1:
                    ri.drlogn[ka-1]+=(gridav[ka])*(log(ni)-log(n0))*iFth/dr
                    ri.drlogn[ka]+=(log(ni)-log(n0))*iFth/dr
                    ri.drlogn[ka+1]+=(1-gridav[ka])*(log(ni)-log(n0))*iFth/dr

                if 1<kb<Nr:
                    ri.drlogn[kb-2]+=-(1-gridav[kb-1])*(log(ni)-log(n0))*iFth/dr
                    ri.drlogn[kb-1]+=-(log(ni)-log(n0))*iFth/dr
                    ri.drlogn[kb]+=-(gridav[kb-1])*(log(ni)-log(n0))*iFth/dr

                ri.dazlogn[ka:kb] += 2*gridav[ka:kb,newaxis]*(log(ni)-log(n0))*dFth

    def write_to_file(self, wl, filename):
        """
        Save annuli to waveguide file readable by Nader's code
        """
        spacer = "-"*80 + "\n"
        import pylab
        numannuli = len(self.annuli)
        minphi = min([a.phi[0] for a in self.annuli])
        maxphi = max([a.phi[1] for a in self.annuli])

        print "Found annuli in range %.6g to %.6g" % (minphi*180/pi, maxphi*180/pi)
        print "Saving %d annuli to '%s'" % (numannuli,filename)

        waveguidefile = open(filename, "w")
        
        waveguidefile.write( spacer )
        waveguidefile.write( "Number of annular sectors\n%d\n" % (numannuli) )
        waveguidefile.write( spacer )
        waveguidefile.write( "real(nhole)\timag(nhole)\tphi[rad]\tdeltaphi[rad]\tr\tdeltar\n")
        waveguidefile.write( spacer )
        
        #Output all annuli to waveguide file
        #Note: Angles are in radians
        for annulus in self.annuli:
            nshape = annulus.material.index(wl)
            ra, rb, phia, phib = annulus.get()

            #pylab.plot([ra*cos(phia),rb*cos(phia),rb*cos(phib)],
            #   [ra*sin(phia),rb*sin(phia),rb*sin(phib)],'k-')
            
            waveguidefile.write( "%.6g\t%.6g\t%.10g\t%.10g\t%.10g\t%.10g\n" \
                % (real(nshape),imag(nshape),phia,(phib-phia),ra,rb-ra) )
        
        waveguidefile.write( spacer+"Number of ellipses\n0\n" )
        waveguidefile.write( spacer+"Resolution of ellipses\n200\n" )
        waveguidefile.write( spacer+"Ellipses data:\n" )
        waveguidefile.write( "real(nhole)\timag(nhole)\tellipticity\trot\txCM\tyCM\tradius\n" )
        waveguidefile.write( spacer+"Radial nodes equally spaced?\nY\n" )
        waveguidefile.write( spacer+"Radial node positions:\n"+spacer)
        
        waveguidefile.close()

# --------------------------------------------------------------------------------------------
#
#  Main waveguide class
#
# --------------------------------------------------------------------------------------------

class Waveguide(object):
    """
    A waveguide object representing the crossectional fiber structure.
    
    Optional arguments are:
        rmin:       minimum radius of the computation (default is currently zero)
        rmax:       maximum radius of the computation (default is automatic)

        material:   a Material object giving the substrate material for the waveguide
        interior:       interior material (defaults to material)
        exterior:       exterior material (defaults to material)
        symmetry:   rotational symmetry of the waveguide
    
        oversample:     how finely sliced the objects are in the azimuthal direction
        reflect:            (Currently not used)
    
    """
    #Proportion of radius to leave between last object and waveguide boundary
    #Defaults the maximum of 2% or 2 computational nodes
    boundary_padding = 0.05
    padding_minimum_nodes = 5
    
    #Don't print information about the creation of the index
    quiet = True
    
    def __init__(self, rmax=None, rmin=0, material=None, interior=None, exterior=None, \
                        symmetry=1, reflect=0, oversample=4):
        self.symmetry = symmetry
        self.reflect = reflect
        self.material = material
        self.oversample = oversample

        #Set interior and exterior materials (defaults to material)     
        self.interior_material = interior
        if exterior is None:
            self.exterior_material = material
        else:
            self.exterior_material = exterior

        #Set rmax manually if specified
        self.set_rrange(rmax, rmin)
        self.set_xyrange()
        
        #Waveguide shape list
        self.shapes = []
    
    def __hash__(self):
        paramhash = int(hash((self.rmax, self.symmetry, self.material)))
        for s in self.shapes:
            paramhash = hash( (paramhash, hash(s)) )
        return paramhash
    
    def __str__(self):
        str = "Waveguide with symmetry=%d, rmax=%.4g, material=%s, " \
                % (self.symmetry, self.rmax, self.material)
        for shape in self.shapes:
            str += "\n> " + shape.__str__()
        return str

    def set_rrange(self, rmax=None, rmin=None):
        "Set manual maximum radius"
        rrange = [None,None]
        if rmax is not None: rrange[1] = rmax
        if rmin is not None: rrange[0] = rmin
        self.rrange_manual = rrange
    def get_rrange(self, boundary_padding=None):
        #Automatic calculation
        rrange = self.extents(boundary_padding)[:2]
        #Check if manual overrides
        if self.rrange_manual[1] is not None: rrange[1] = self.rrange_manual[1]
        if self.rrange_manual[0] is not None: rrange[0] = self.rrange_manual[0]
        return rrange

    def get_xyrange(self):
        "Return the boundary domain size for cartesian calculations"
        if self.xyrange_manual is not None:
            return self.xyrange_manual
        else:
            return self.extents(xy=True)
    def set_xyrange(self, xyrange=None):
        "Set manual x and y limits"
        self.xyrange_manual = xyrange

    def set_rmax(self, rmax, boundary_padding=None):
        self.set_rrange(rmax)
    def get_rmax(self, boundary_padding=None):
        return self.get_rrange(boundary_padding=boundary_padding)[1]
    rmax = property(get_rmax, set_rmax)
    
    def calc_core_size(self):
        "Guess approximate core size from shapes in waveguide"
        if hasattr(self, 'coresize_manual'):
            return self.coresize_manual
        
        #Guess the smallest radius in the fiber
        rc = self.extents(0)[0]
        
        #If a shape encloses the origin find minimum maximal radii
        if rc<=0:
            rc = inf
            for s in self.shapes:
                rmin, rmax, phimin, phimax = s.extents()
                if rmin<=0: rc = min(rc, rmax)
        return rc
    def set_core_size(self, rc):
        self.coresize_manual = rc
    core_size = property(calc_core_size, set_core_size)

    #Routines to calculate and set the interior and exterior  material
    def _calc_interior_material(self):
        "Guess the material of the core from all shapes in waveguide"
        if self.interior_material_manual is not None:
            return self.interior_material_manual
        
        #Default to waveguide material
        core_material = self.material
        
        #Find any shapes enclosing the origin
        #or shapes touching minimum radius
        for s in self.shapes:
            rmin, rmax, phimin, phimax = s.extents()
            
            if rmin<=0 and not isinstance(s, Image):
                core_material = s.material

        return core_material
    def _set_interior_material(self, mat):
        self.interior_material_manual = mat
    interior_material = property(_calc_interior_material, _set_interior_material)
    
    def extents(self, boundary_padding=None, xy=False):
        "Return the maximum radius in the numerical solution"
        #Manual override for cartesian domain
        if not xy and self.xyrange_manual is not None:
            return self.xytange_manual

        bmin = [inf]*4; bmax=[-inf]*4
        for shape_ in self.shapes:
            bounds = shape_.extents(xy=xy)
            bmin = minimum(bmin,bounds)
            bmax = maximum(bmax,bounds)
        
        #Use internal padding
        if boundary_padding is None:
            boundary_padding = self.boundary_padding

        #Add padding to inner and exterior limits           
        Nround = 256
        if xy and boundary_padding is not None:
            offsetx = ceil(bmax[1]*boundary_padding*Nround)/Nround
            offsety = ceil(bmax[3]*boundary_padding*Nround)/Nround
            bmax[1] += offsetx; bmin[0] -= offsetx
            bmax[3] += offsety; bmin[2] -= offsety
        elif not xy and boundary_padding is not None:
            offsetr = ceil(bmax[1]*boundary_padding*Nround)/Nround
            
            #Don't add if external material is defined or different from material
            if self.exterior_material is self.material:
                bmax[1] += offsetr
            if self.interior_material is self.material:
                bmin[0] = max(bmin[0]-offsetr,0)

        #Manual overrides for radial domain                 
        if not xy and self.rrange_manual[0] is not None:
            bmin[0] = self.rrange_manual[0]
        if not xy and self.rrange_manual[1] is not None:
            bmax[1] = self.rrange_manual[1]

        return [ bmin[0], bmax[1], bmin[2], bmax[3] ]

    def get_coord(self, Nshape, border=0):
        rmin, rmax = self.get_rrange()

        # Use a pretty high value for the calculation bandwidth to give good
        # accuracy in the mode calculations
        dbw = 7
        
        return WgCoord(rmax, rmin, symmetry=self.symmetry,
                reflect=self.reflect, Nshape=Nshape, border=border, bandwidth=dbw)

    # +-----------------------------------------------------------------------+
    # | Pickling functions
    # | We don't pickle the cached data to give compact storage
    # +-----------------------------------------------------------------------+
    def __getstate__(self):
        "Pickle all needed data, ignore cached data"
        state = self.__dict__.copy()
        #Remove cached information
        for key in state.keys():
            if key.startswith('_'):
                del state[key]
        return state
    
    def __setstate__(self,state):
        "Restore needed data"
        self.__dict__.update(state)
    
    # +-----------------------------------------------------------------------+
    # | Shape and Waveguide functions
    # +-----------------------------------------------------------------------+
    def add_shape(self, *shapes):
        "Add a shape or shapes to the waveguide"
        for shape in shapes:
            #If we are given an array of shapes, add them individually
            if iterable(shape):
                self.add_shape(*shape)
                continue

            #Check that it is a shape before adding
            if not isinstance(shape, WgShape):
                logging.error( "%s doesn't seem to be a Waveguide Shape!" % (shape) )
    
            self.shapes.append( shape )

    add_shapes = add_shape

    def remove_shape(self, *shapes):
        "Remove a shape or shapes from the waveguide"
        if len(shapes)==0:
            self.shapes = []
        else:
            for shape in shapes:
                self.shapes.remove(shape)
        
    def calc_paramhash(self, *params):
        "Calculate hash for waveguide, including shapes and parameters"
        paramhash = [ hash(self.rmax), hash(self.symmetry), \
            hash(self.reflect), hash(self.material) ]
        paramhash += [ hash(p) for p in self.shapes ]
        paramhash += [ hash(p) for p in params ]
        return paramhash

    def clear(self):
        if hasattr(self, '_combined_'):
            del self._combined_

    # +-----------------------------------------------------------------------+
    # | Index functions
    # +-----------------------------------------------------------------------+
    def calculate(self, coord, wl, mask=None, resample=None):
        '''
        Calculates the refractive index for the grid of shape Nshape iterating
        over all included shapes. Returns a combined index object.
        
        In most cases use the helper methods: index, index2, dlognr and dlognaz
        
        Nshape: (Nr,Naz) to use
        wl: wavelength to calculate refractive indicies at
        border: n to eliminate n border points at either end
        mask: return a mask which is 1 where the material is the same as that give in the mask
        '''
        #If we have already calclated the index with these parameters, return the
        #cached calculation
        paramhash = self.calc_paramhash(coord, wl, mask)
        if hasattr(self, '_combined_') and self._combined_hash_==paramhash:
            return self._combined_

        if not self.quiet:
            logging.info( "Generating Waveguide from %d object%s; rmax:%.2g, Nr:%d, Naz:%d" \
                % (misc.numands(len(self.shapes)) + (self.rmax, coord.Nr, coord.Naz)) )

        #Check boundary padding
        rshapes = self.get_rmax(boundary_padding=0)
        rmax = self.get_rmax()
        nodes_in_padding = floor((rmax-rshapes)/coord.dr)
        if not self.quiet and nodes_in_padding<self.padding_minimum_nodes:
                logging.warning( "Padding nodes %d less than minimum, increase rmax?" \
                % nodes_in_padding )

        #Combined object to collect n(r,φ)² d/dr log[n(r,φ)²] and d/dφ log[n(r,φ)²]
        if mask is not None:
            combinedri = CombinedMask(coord, self.material, mask, dtype=complex_)
        else:
            combinedri = CombinedIndex(coord, self.material, wl, dtype=complex_)

        #Combine the shapes into a single discretized refractive index. Use various
        # types of merging based on what methods the shapes expose,
        # currently either annular or index based methods.
        amerge = AnnularCombine(coord, self.symmetry)
        dmerge = DirectCombine(coord)
        for shape in self.shapes:
            if hasattr(shape, "get_annuli"):
                amerge.merge_shapes([shape], oversample = self.oversample)
        
            elif hasattr(shape, "get_index"):
                #Does not yet support oversampling
                dmerge.merge_shapes([shape])
        
        amerge.generate(combinedri)
        dmerge.generate(combinedri)
        
        #Cache the combine object to speed sucessive accesses (with the same paramters)
        self._combined_ = combinedri
        self._combined_hash_ = paramhash
        return combinedri

    # +-----------------------------------------------------------------------+
    # | Index functions
    # +-----------------------------------------------------------------------+
    def calculate_radial_layers(self, wl, mask=None):
        '''
        Calculates the refractive indices and radii assuming a azimuthally independent
        shapes.
        '''
        if not self.quiet:
            logging.info("Generating index profile from %d object%s" % misc.numands(len(self.shapes)))
        
        from LayeredSolver import MidLayer, InteriorLayer, ExteriorLayer
        
        #Add layers from shapes
        layer_list = []
        for shape in self.shapes:
            r1,r2 = shape.extents()[:2]
            n = shape.material.index(wl)
            z = shape.zorder

            layer_list.append(MidLayer(r1,r2,n,zorder=z))
        
        #Sort on rstart for each layer
        layer_list.sort()
        
        #Process the list and create true layers
        nint = self.interior_material.index(wl)
        n0 = self.material.index(wl)
        next = self.exterior_material.index(wl)

        new_layers=[]
        last=None
        for layer in layer_list:
            #Add first layer from wg or shape
            if last is None:
                if layer.r1>0:
                    new_layers.append(InteriorLayer(0,layer.r1,nint))
                else:
                    layer = InteriorLayer(0,layer.r2,layer.n,zorder=layer.zorder)

            #Add internal layer
            else:
                #If there is a gap insert a new layer
                if layer.r1>last.r2:
                    new_layers.append(MidLayer(last.r2,layer.r1,n0))
                
                #Fix overlapping annulii
                elif layer.r1<last.r2 and layer.r2>=last.r2:
                    #Fix the layer transition appropriately
                    if layer.zorder>=last.zorder:
                        last.r2=layer.r1
                    else:
                        layer.r1=last.r2

                elif layer.r1<last.r2 and layer.r2<last.r2:
                    #Fix the layer transition appropriately
                    if layer.zorder>=last.zorder:
                        last.r2=layer.r1                            #Fix last layer
                        new_layers.append(layer)        #Add this layer
                        layer = MidLayer(layer.r2,last.r2,last.n,zorder=last.zorder)
                    else:
                        layer=None
            
            #Add this layer if it still exists
            if layer is not None:
                new_layers.append(layer)
                last=layer

        #Final exterior layer - also check for chirality!
        if last.n==next:
            new_layers[-1] = ExteriorLayer(last.r1,last.r2,last.n,zorder=last.zorder)
        else:
            new_layers.append(ExteriorLayer(layer.r2,inf,next))
        return new_layers

    def mask(self, Nshape=None, mat=None, coord=None, resample=None, border=0):
        '''
        Calculates a material mask returning 1 wherever the WG contains the given
        material and 0 otherwise.
        '''
        if coord is None: coord = self.get_coord(Nshape, border=0)
        
        if mat is None:
            mat = self.material
        
        #Calculate mask object
        nix = self.calculate(coord, 1, mask=mat)
        
        #Resample to another coordinate
        if resample:
            ext = (self.interior_material==mat, self.exterior_material==mat)
            m2 = resample.fourier_resample(fftp.fft(nix.ni2, axis=1), coord, ext=ext)
        else:
            m2 = nix.ni2[self.borderslice(border)]
        return m2

    def exterior_index(self, wl):
        "The exterior (r>rmax) substrate refractive index"
        return self.exterior_material.index(wl)
    index0 = exterior_index

    def interior_index(self, wl):
        "The interior (r<rmin) substrate refractive index"
        return self.interior_material.index(wl)

    def index_range(self, wl):
        "The minimum and maximum refractive index in the wg"
        minindex=self.material.index(wl)
        maxindex=self.material.index(wl)
        for shape in self.shapes:
            sindex = shape.material.index(wl)
            minindex = min(minindex, sindex)
            maxindex = max(maxindex, sindex)
        return minindex,maxindex

    def borderslice(self,border):
        if border>0:
            bslice = slice(border, -border)
        else:
            bslice = slice(None)        
        return bslice

    #Convenience functions:
    def index(self, *args, **kwargs):
        ni = sqrt(self.index2(*args, **kwargs))
        return ni

    def index2(self, wl, Nshape=None, coord=None, resample=None, border=0):
        if coord is None: coord = self.get_coord(Nshape, border=0)
        
        #Calculate index object
        nix = self.calculate(coord, wl)
        
        #Resample to another coordinate
        if resample:
            ext = (self.interior_material.index(wl)**2, self.exterior_material.index(wl)**2)
            ni2 = resample.fourier_resample(fftp.fft(nix.ni2, axis=1), coord, ext=ext)
        else:
            ni2 = nix.ni2[self.borderslice(border)]
        return ni2

    # +-----------------------------------------------------------------------+
    # | Plotting and output methods
    # +-----------------------------------------------------------------------+

    ## Plot the waveguide .. must provide resolution
    def plot_bitmap(self, Nx=(200,21), coord=None, average=False, sectors=None, wl=1, **opts):
        '''
        Plots the refractive index of the waveguide and included shapes at the
        given wavelength. The wavelength defaults to 1.

        Supply one of:
        Nx:             number of radial & azimuthal nodes
        coord:      a coordinate object on which the refractive index is sampled
        
        Optionally:
        average: azimuthal averaging of refractive index
        sectors: number of symmetry sectors to plot
        wl:         wavelength to evaluate the refractive index
        
        Other options:
        Other Plotting style options can be given
        style:          pcolor, 3d, line and contour
        cmap:       color map from matplotlib.cm
        '''
        from Plotter import plot_v, title, colorbar, pcolor, plot, xlabel, ylabel, legend, draw

        if opts=={}: opts = {'style':'pcolor'}
        if coord is None: coord = self.get_coord(Nx, border=0)

        #Calculate the refractive index
        inx = fft.fftshift(self.index(wl, coord=coord), axes=[1])

        #Refractive index of all materials in waveguide
        materials = [shape.material for shape in self.shapes] + [self.material]
        ticks = [real(mx.index(wl)) for mx in unique(materials)]
        
        #Average the refractive index azimuthally
        if average: inx = inx.mean(axis=1)
        
        #Handle 1d-plotting
        if inx.ndim==1 or inx.shape[1]==1:
            plot(coord.rv, inx.real, '-', label="Real")
            
            #Only plot imaginary part if large enough
            if absolute(imag(inx)).max()>1e-8:
                plot(coord.rv, inx.imag, '--', label="Imag")
                
            #Plot material indicies
            for tick in ticks:
                plot([coord.rv.min(), coord.rv.max()],[tick,tick],':', label='%.5g'%real(tick))
            
            ylabel("Refractive index")
            xlabel('Radial distance, r')
            legend(loc='best')
        else:
            #Plot all waveguide sectors by default
            if sectors is None: sectors = self.symmetry
            coord.sectors = sectors

            #Paint in multiple sectors simply
            final_inx = inx.copy()
            for sec in range(sectors-1):
                final_inx = append(final_inx, inx, axis=1)

            #Plot the data!
            try:
                plot_v(coord, real(final_inx), **opts)
        
                title("Refractive index")
                colorbar(fraction=0.05, pad=0.025, spacing='uniform', format='%.4f', ticks=ticks)
            except NotImplementedError:
                logging.error("Plot style not recognised or plotting failed")

        draw()
        
    def plot(self, rotate=0, substrate=True, fill=True, sectors=None, **user_style):
        '''
        Plot the shapes in the waveguide
        
        Optional parameters:
         rotate:            rotate plot by specified radians
         substrate:     plot waveguide boundary
         sectors:           number of sectors
         fill:                  Plot the shapes filled or as outline, True or False
         colour:            Plot all objects with one color
        
        Other arguments:
        Any matplotlib style commands, for example:
        * fc = foreground color
        * ec = edge color
        * lw = linewidth
        '''
        import pylab as p_
        from matplotlib import patches

        #Number of sectors to plot
        if sectors is None: sectors = self.symmetry
        
        if sectors<self.symmetry: substrate = False

        #Default style
        style = {}

        #Override edge and fill colors
        if 'colour' in user_style:
            user_style['ec'] = user_style['fc'] = user_style.pop('colour')
        if 'color' in user_style:
            user_style['ec'] = user_style['fc'] = user_style.pop('color')

        logging.info( "Plotting Waveguide from %d object%s" \
                % misc.numands(len(self.shapes)) )

        ax = p_.gca()

        #Plot the boundary first!
        if substrate:
            style={'fc':self.material.color, 'ec':'black', 'ls':'dashed', 'fill':fill, 'alpha':0.5, 'lw':(not fill)*1}
            style.update(user_style)
            
            rmin,rmax = self.get_rrange()
            ax.add_artist(patches.Circle((0,0), radius=rmax, **style))
            ax.add_artist(patches.Circle((0,0), radius=rmin, **style))
    
        #Draw symmetry objects for each shape
        for shape in sorted(self.shapes):
            logging.debug( "Plotting %s at %s" % ((shape.name, shape.center)) )

            #Color shape based on material, unless overrided
            style = {'ec':shape.material.color, 'fc':shape.material.color, 'fill':fill, 'lw':(not fill)*1}
            style.update(user_style)

            #Plot the shape for each of the sectors
            for nx in range(0,sectors):
                theta=nx*2*pi/self.symmetry + rotate
                objects = shape.get_plot_object(theta, **style)
                
                #Add all objects in plot shape
                for o in objects: ax.add_artist(o)
    
        #We need to control the axis ourselves
        ax.set_xlabel(r'$x,\ \mu$m')
        ax.set_ylabel(r'$y,\ \mu$m')
        ax.axis([-self.rmax,self.rmax,-self.rmax,self.rmax])
        ax.set_aspect('equal')

        p_.draw()

