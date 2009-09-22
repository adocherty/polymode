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
Image loading & regridding

To Do:
* Add other image loaders .. currently requires matplotlib.imread
* 
"""

from __future__ import division

import logging
from numpy import *

def rms(f,x):
    rms = trapz(trapz(x*absolute(f)**2)) / trapz(trapz(absolute(f)**2))
    return rms

def bracket(x):
    return minimum(1.0,maximum(0.0,x))

#Class for array with index checking making it easier to perform interpolation
class guardedarray(ndarray):
    def __getitem__(self,index):
        try:
            out = ndarray.__getitem__(self, index)
        except IndexError:
            if hasattr(index, '__len__'):
                ilen = len(index)
            else:
                ilen = 1
                index=[index]
            
            if self.ndim<>ilen:
                raise IndexError, "Index has more dimentions than array"

            for dim in xrange(ilen):
                if isinstance(index[dim], ndarray):
                    index[dim][index[dim]>=self.shape[dim]] = self.shape[dim]-1
                elif isinstance(index[dim], slice):
                    pass
                elif isinstance(index[dim], int):
                    index[dim] = max(min(index[dim],self.shape[dim]),0)
                    
            out = ndarray.__getitem__(self, index)
        return out

def create_grid(r_max,theta_max,Nr,Ntheta):
    dr = r_max/Nr; dtheta = theta_max/Ntheta
    rv = arange(0,Nr+2)*dr
    thetav = arange(0,Ntheta)*dtheta
    kv = fft.fftshift(arange(-Ntheta//2,Ntheta//2))*2*pi/theta_max
    return dr,dtheta,rv,thetav,kv

#Regrid using nearest neighbour interpolation
def nn_polar_regrid(density,center,r,theta, data):
    xpolar = r[:,newaxis]*cos(theta)
    ypolar = r[:,newaxis]*sin(theta)
    xarray,yarray = floor(density[0]*xpolar+center[0]).astype(int),floor(density[1]*ypolar+center[1]).astype(int)

    polar_data = data[xarray, yarray]
    return polar_data

#Regrid using bilinear filtering
def bilinear_polar_regrid(density,center,r,theta, data):
    xpolar = density[0]*r[:,newaxis]*cos(theta)+center[0]
    ypolar = density[1]*r[:,newaxis]*sin(theta)+center[1]
    xarray,yarray = floor(xpolar).astype(int), floor(ypolar).astype(int)
    xfrac,yfrac = xpolar%1, ypolar%1
    
    polar_data = data[xarray, yarray]*xfrac*yfrac + data[xarray+1, yarray]*(1-xfrac)*yfrac \
        + data[xarray, yarray+1]*xfrac*(1-yfrac) + data[xarray+1, yarray+1]*(1-xfrac)*(1-yfrac)
    return polar_data

#Regrid using inverted nearest neighbour binning
def nn_inverse_polar_regrid(density,center,r,theta, data):
    xpolar = (arange(data.shape[0])-center[0])/density[0]
    ypolar = (arange(data.shape[1])-center[1])/density[1]

    rx = sqrt(xpolar[:,newaxis]**2 + ypolar[newaxis,:]**2)
    thetax = arctan2(ypolar[newaxis,:], xpolar[:,newaxis])
    
    dr = r[1]-r[0]; dtheta = theta[1]-theta[0]
    rmin,rmax = r[0],r[-1]; tmin,tmax = theta.min(), theta.max()

    polar_data = zeros((len(r), len(theta)), float_)
    count = zeros((len(r), len(theta)), float_)

    #iterate over all pixels in source image and map them to the polar image
    #by nearest neighbour
    for ii in range(data.shape[0]):
        for jj in range(data.shape[1]):
            if rx[ii,jj]<rmax and tmin<thetax[ii,jj]<tmax:
                pinx = (floor((rx[ii,jj]-rmin)/dr+0.5),floor((thetax[ii,jj]-tmin)/dtheta+0.5))
                polar_data[pinx] += data[ii,jj]
                count[pinx] += 1
    
    #If count==0 then take an average of the surrounding pixels
    for ii,jj in transpose(nonzero(count==0)):
        xp = density[0]*r[ii]*cos(theta[jj])+center[0]
        yp = density[1]*r[ii]*sin(theta[jj])+center[1]

        xa,ya = floor(xp).astype(int), floor(yp).astype(int)
        xf,yf = xp%1, yp%1

        polar_data[ii,jj] = data[xa, ya]*xf*yf + data[xa+1, ya]*(1-xf)*yf \
            + data[xa, ya+1]*xf*(1-yf) + data[xa+1, ya+1]*(1-xf)*(1-yf)
        count[ii,jj] = 1

    polar_data/=count
    return polar_data


#Calculate the overlap difference in the polar image with symmeties cs
def calc_symmetry_factor(f_polar,kv,cs):
    Nr=shape(f_polar)[0]

    #attach more weight to the information from smaller r
    r_weight = 1 #exp(-arange(Nr)[:,newaxis]/Nr)
    f_fourier_sum = fft.fft( average(absolute(f_polar*r_weight), axis=0) )
    if isinstance(cs,ndarray):
        return sum( absolute(f_fourier_sum) * sin(pi*kv/cs[:,newaxis])**2, axis=1)/sum(sin(pi*kv/cs[:,newaxis])**2)
    else:
        return average( absolute(f_fourier_sum), weights=sin(pi*kv/cs)**2 )

#Calculate the fitness of the polar regridding with a particular symmetry and center
def find_symm_factor(center, f, coord, density, polar_regrid=bilinear_polar_regrid):
    rv, phiv = coord.rv, coord.calc_phiv(sectors=coord.symmetry)
    kv = coord.calc_ms(sectors=coord.symmetry)
    f_polar = polar_regrid(density, center, rv, phiv, f)
    return calc_symmetry_factor(f_polar,kv,coord.symmetry)


class ImportBitmap:
    '''
    Import a bitmap image (png) file

    Parameters:
     - filename: name of image file (not file object)
     - density: dpi of image either single numer or tuple of (dpi_x, dpi_y)
     - size: size in um of image either single numer or tuple of (size_x, size_y)
     this is used to calculate the density, so only one of density, size should be given
     - hist: histogram function to use to normalize the image data
     - negative
        if False black is interpreted as no material, white is interpreted as presense of material
        if True the opposite interpretation is used
    '''
    auto_adjust=True
    histogram_levels=64
    def __init__(self, filename, density=None, size=None, rmax=None, center=None, offset=0, negative=False):

        try:            #Try using Python Imaging Libraries to load image
                    from PIL import Image
                    self.image = Image.open(filename)
                    idata = array(self.image).astype(float32)
        except:         #Fallback to matplotlib libraries (soon to be LodePNG?)
                    from pylab import imread
                    idata = imread(filename)
                    
        #If RGB then average over all components
        #This isn't of course the perceptual intensity
        if idata.ndim>2:
            idata = idata.mean(axis=2)

        #Form a guarded array to prevent exceptions when accessing out of area
        self.image_data = guardedarray(idata.shape, buffer=idata, dtype=idata.dtype)
        
        #Invert the image
        if negative: self.image_data = -self.image_data

        #Set to range between 0 and 1
        self.normalize(self.image_data, clip=False, hist=False)
        
        #Other information
        self.Nx,self.Ny = shape(self.image_data)
        self.center = center
        self.offset = offset
        self.polar_regrid=nn_polar_regrid

        if not size is None:
            density = self.Nx/size[0], self.Ny/size[1]
        if density is None:
            density = 10.,10.
        self.density = density
        if rmax is None:
            self.rmax = min(self.Nx/density[0], self.Ny/density[1])/2
            logging.debug("Auto rmax is %.2g" % self.rmax)
        else:
            self.rmax = rmax
        
        logging.info( "Loaded image file %s" % filename )

    def get_center(self, coord=None):
        '''
        Find the optimal center that gives the best symmetry of the image
        shape: (Nr,Naz) the number of points needed in the final image
        symmetry: symmetry of image
        rmax: only use area within r<rmax of the center
        '''
        approx_center = self.Nx/2., self.Ny/2.
        
        if self.center is not None:
            center = self.center
        elif coord is None or coord.symmetry==1:
            center = approx_center
        else:
            from scipy import optimize

            symmetry = coord.symmetry
            shape = coord.Nr, coord.Naz
            rmax = coord.rmax

            center=optimize.fmin( find_symm_factor, approx_center, 
                args=(self.image_data, coord, self.density,self.polar_regrid), disp=0)
            logging.info("Found optimal center: x=%.6g, y=%.6g" % tuple(center))
        
        return center
        
    def regrid(self, coord, sectors=None, rborder=0):
        '''
        Perform cartesian to polar regridding on the loaded image
        sectors: average over this many sectors (1<sectors<symmetry)
        rborder: add a border of this many rows to the final radial data of constant value
        '''
        if sectors is None:
            sectors = coord.symmetry
        center = self.get_center(coord)

        logging.info( "Regridding image to polar coordinates from %d sectors" % sectors )

        #dr,dtheta,rv,thetav,kvx = create_grid(rmax,(2*pi/symmetry)*sectors,Nr,Ntheta*sectors)
        Nr, Naz = coord.Nr, coord.Naz
        phiv = coord.calc_phiv(sectors=sectors)
        polar_data = self.polar_regrid(self.density, center, coord.rv, \
                        phiv+self.offset, self.image_data)

        polar_data = polar_data.reshape( (Nr+2-2*coord.border, sectors, Naz) ).mean(axis=1)

        #Second normalize without histogram to insure max index is 1
        if self.auto_adjust:
            self.normalize(polar_data, clip=True, hist=True)
        else:
            self.normalize(polar_data, hist=False)
        
        # Add a border:
        if rborder>0:
            polar_data[-rborder:] = 0
        
        return polar_data

    def contrast(self, scale=1, transition=0.5, clip=False):
        '''
        scale: The strength of the contrast stretching
        transition: From 0 to 1, the center of the strectch
        '''
        imin = self.image_data.min()
        imax = self.image_data.max()
        tpoint = transition*(imax-imin)+imin
        
        self.image_data = 0.5*tanh(scale*(self.image_data-tpoint))+0.5
        self.image_data[self.image_data<0.01*imin] = imin
        self.image_data[self.image_data>0.99*imax] = imax

    def normalize(self, im, clip=False, hist=False):
        if hist:
            dist,value = histogram(im, bins = self.histogram_levels)
            dist_sort = dist.argsort()
            white = value[ dist_sort[-10:] ].max()
            black = value[ dist_sort[-10:] ].min()
        else:
            black = im.ravel().min()
            white = im.ravel().max()
    
        logging.debug("Normalizing: black=%.2f, white=%.2f" %(black,white))
        im -= black
        im /= (white-black)

        if clip and hist:
            im[:] = bracket( im )
        return im

    def gamma(self, gamma=0.5):
        "Change the gamma of the image."
        imin = self.image_data.min()
        imax = self.image_data.max()
        self.image_data = (self.image_data-imin)**gamma*(imax-imin)**(1-gamma)+imin

    def plot_histogram(self):
        import pylab
        dist,value = histogram(self.image_data, self.histogram_levels)
        dist_sort = dist.argsort()
        white = value[ dist_sort[-10:] ].max()
        black = value[ dist_sort[-10:] ].min()

        pylab.plot(value, dist/dist.max(), 'b:')
        pylab.plot([black,black],[0,1],'k--')
        pylab.plot([white,white],[0,1],'r--')
        
    def plot(self):
        from pylab import clf, imshow, cm, plot
        
        center = asarray(self.get_center())
        psize = asarray(self.image_data.shape, dtype=float)

        lower = -center/asarray(self.density)
        upper = (psize-center)/asarray(self.density)

        clf()
        imshow(self.image_data,extent=(lower[0],upper[0],lower[1],upper[1]), cmap=cm.gray)
        plot([0], [0], 'r+')


