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
Materials for use in constructing Waveguide shapes.
Each material has a refractive index that is dependent upon
the wavelength and (in future) the temperature

 ToDo:
 *  Get info from data files
 *  Restrict plot to wavelength limits for interpolated materials
 *  Add temperature dependence
 *  Use limits for all materials, not just interpolated ones
"""

from __future__ import division
import logging

from numpy import (pi, asarray, array, arange, isscalar, vectorize, real, imag,
    newaxis, any, size, sum, squeeze, shape, ones, iterable, absolute, argsort, atleast_1d,
    nonzero, zeros)
from numpy.lib.scimath import *

from .difflounge import finitedifference
from .mathlink import constants, utf8out

aa_=asarray
material_resource_dirs = [ 'Sopra', 'freesnell', 'luxpop' ]

class Material:
    color = '#656565'
    auto_vectorize = True
    wl_limits = (0.1,10)        #in microns
    temp_limits = (0,400)   #in Kelvin
    T = 300
    citation = ""

    def __str__(self):
        if hasattr(self, 'name'):
            try:
                dscr = self.name.encode('utf-8')
            except:
                dscr = self.__class__.__name__
        else:
            dscr = self.__class__.__name__
        return dscr
        
    def info(self):
        print "%s" % self
        print "Wavelength limits (um): %.4g->%.4g" % tuple(self.wl_limits)
        print "Temperature limits (K): %.4g->%.4g" % tuple(self.temp_limits)
        print "Cite: %s" % self.citation
        
    def wavelength_derivative(self, x, dx=0.01, units='wl'):
        '''
        Numerically evaluate the derivative of the refractive index squared
        for this material at a specific wavenumber.
        '''
        #Evaluate by 5 point FD around wavelength
        Nb = 5
        wavelength = self.convertwavelength(x, units)
        wls = arange(-(Nb//2),(Nb+1)//2)*dx + wavelength
        n2 = self.index(wls, 'wl')**2
        
        Dl = finitedifference.FiniteDifference(bandwidth=Nb, X=wls)
        return Dl.diff_at(n2, wavelength)

    def index(self, x, units='wl'):
        '''
        Get the refractive index for this material at a specific wavenumber.
        '''
        wavelength = self.convertwavelength(x, units)
        
        if not hasattr(self, 'limits_warning_shown'):
            if any(wavelength<min(self.wl_limits)) or any(wavelength>max(self.wl_limits)):
                if size(wavelength)>1:
                    wlrange = "%.4g-%.4g" % (min(wavelength), max(wavelength))
                else:
                    wlrange = "%.4g" % (wavelength)

                logging.warning( utf8out(u"The wavelength %sμm is outside the valid\n" % (wlrange)
                + u"\trange for this material, large errors may result.\n"
                + u"\t%s only has data in the range (%.4g, %.4g)μm" % ((self,)+self.wl_limits)) )
                
#           if any(self.T<min(self.temp_limits)) or any(self.T>max(self.temp_limits)):
#               logging.warning("Refractive index is being requested outside the valid temperature")
#               logging.warning("range for this material, large errors may result.")
#               logging.warning("%s only has data in the range (%.4g, %.4g)K" 
#                   % ((self,) + self.temp_limits) )
            self.limits_warning_shown = True
            
        if iterable(wavelength) and self.auto_vectorize:
            ifunction = vectorize(self.index_function)
        else:
            ifunction = self.index_function
        ans = ifunction(wavelength)
        return ans

    def convertwavelength(self, x, wlunits='um'):
        '''
        Convert wavelength/wavenumber to wavelength in um
        '''
        if wlunits=='ev':                           #electronvolts
                                                wl = constants.h*constants.c/constants.eV*1e6/x
        elif wlunits=='/cm':                    #inverse cm
                                                wl = (1e-2/x) * 1e6
        elif wlunits=='nm':                         #nanometers
                                                wl = 1e-3*x
        elif wlunits=='wn':                         #Wavenumber
                                                wl = 2*pi/x
        elif wlunits=='A':                          #Angstroms
                                                wl = 1e-4*x
        elif wlunits=='um' or wlunits=='wl':        #micrometers
                                                wl = x
        else:
            wl = x
        
        return wl

    def plot(self,wlrange=None, points=100, showdata=False):
        '''
        Plot the material refractive index
        * wrange is the wavelength range
        * points is the number of points to use in the plot
        '''
        from . import Plotter as pl
        
        #Plot over supplied range, internal wl limits or
        #an arbitrary default
        if wlrange is None: wlrange = self.wl_limits

        dwl = float(wlrange[1]-wlrange[0])/points
        wls = arange(points)*dwl + wlrange[0]+dwl/2
        nis = self.index(wls)
        
        #Plot real part
        pl.plot(wls, real(nis), '-', label="%s [re]" % self)
        
        #Only plot imaginary part if large enough
        if absolute(imag(nis)).max()>1e-8:
            pl.plot(wls, imag(nis), '--', label="%s [im]" % self)

        pl.ylabel(r"$\rm{Refractive Index},\ n$")
        pl.xlabel(r"$\rm{Wavelength},\ \mu m$")
        pl.legend(loc='best')
        pl.grid(True)

#----------------------------------- Base Material Classes  ------------------------------------#

class Fixed(Material):
    def __init__(self,ni):
        self.ni = ni
    def index_function(self, wavelength):
        return self.ni
    def __str__(self):
        return "Fixed Index (n=%s)" % self.ni

class CauchyApproximation(Material):
    '''
    CauchyApproximation(A=[])
    A material with refractive index calulated from the Cauchy approximate given by
    I. Nikolov and C. Ivanov in Proc. SPIE, Vol. 3573, 409 (1998); doi:10.1117/12.324553
    The coefficients should be supplied in array A.
    
        ni^2 = A(0) + A(1) l^2 + A(2) l^(-2) + A(3) l^(-4) + ...
    '''
    A=[]
    def __init__(self, A=None):
        if A is not None:
            self.A = A

    def index_function(self, wavelength):
        Lterms = len(self.A)
        lpows = -2*arange(-1,Lterms-1)
        
        #First two indices are swapped!
        lpows[:2]=[0,2]
        
        return sqrt(sum(self.A*wavelength**lpows))

class Sellmeier(Material):
    '''
    Sellmeier(B=[], L=[], A=1)
    A material with refractive index calulated from the Sellmeier equation
    The coefficients should be supplied in arrays B and L. Each pair of coefficients represents
    an absorption resonance of strength Bi at a wavelength Li (in µm).
    
    The Sellmeier equation:
        ni^2 = A + B(0) l^2 / (l^2 - L(0)^2 + B(1) l^2 / (l^2 - L(1)^2 + ..
    '''
    A = 1
    def __init__(self, B=None, L=None, A=None):
        if B is not None and L is not None:
            self.B = B
            self.L = L
    def index_function(self, wavelength):
        eps = self.A + sum( aa_(self.B)*wavelength**2 / (wavelength**2-aa_(self.L)**2))
        return sqrt(squeeze(eps))

class ClassiusMorlotti(Material):
    '''
    ClassiusMorlotti(B=[], L=[])
    A material with refractive index calulated from the Classius-Morlotti equation
    [citation needed]
    
    The Classius-Morlotti equation:
        (n^2-2)/(n^2+2) = B(0) l^2 / (l^2 - L(0)^2 + B(1) l^2 / (l^2 - L(1)^2 + ..
    '''
    def __init__(self, B=None, L=None):
        if B is not None and L is not None:
            self.B = B
            self.L = L

    def index_function(self, wavelength):
        if hasattr(self, 'D'):
            eps = sum( aa_(self.B + self.D*self.f)*wavelength**2
                    / (wavelength**2-aa_(self.L)**2))
        else:
            eps = sum( aa_(self.B)*wavelength**2 / (wavelength**2-aa_(self.L)**2))
        
        return squeeze(sqrt((2*eps+1)/(1-eps)))

class Laurent(Material):
    '''
    A material with refractive index calulated from the Laurent series
    The coefficients of the series should be supplied in array A. The wavelength
    is in µm.
    '''
    def __init__(self, A=None):
        if A is not None:
            self.A = A
    def index_function(self, wavelength):
        series_l = wavelength**(2*(1-arange(0,len(self.A))))
        eps = sum( aa_(self.A)*series_l )
        return sqrt(squeeze(eps))

class Composite(Material):
    '''
    Composite(materials = [], weights = [])
    A material with refractive index calulated from the weighted average of two or more 
    previously defined materials.
    '''
    def __init__(self, materials = [], weights = [], interpolation='CM'):
        self.materials = materials
        self.weights = weights
        
        self.wl_limits = (max([ min(m.wl_limits) for m in materials ]),
                        min([ max(m.wl_limits) for m in materials ]))
                        
        self.itype = interpolation
    
    def info(self):
        print "Composite of:"
        for mat,w in zip(self.materials, self.weights):
            print "> %.2g%% \t%s " % (w*100, mat)
        
    def index_function(self, wavelength):
        Nmaterial = len(self.materials)
        #Weighted or regular average
        if len(self.weights)==Nmaterial:
            weights = self.weights
        else:
            weights = [1./Nmaterial,]*Nmaterial
        
        #Average the materials using linear or Classius-Morlotti type interpolation
        ni = 0
        for ii in range(Nmaterial):
            if self.itype=='linear':
                ni += self.materials[ii].index_function(wavelength)**2 * weights[ii]
            else:
                n = self.materials[ii].index_function(wavelength)
                ni +=  (n**2-1)/(n**2+2) * weights[ii]
        
        if self.itype=='linear':
            ni = sqrt(ni)
        else:
            ni = sqrt((2*ni+1)/(1-ni))
        
        return ni

class Interpolated(Material):
    '''
    A material with refractive index calulated by interpolating from a list
    of refractive indices at different wavelengths.
    Although there are wavelength limits which warn if you are outside the
    range of data, some datasets can have large gaps in which the values
    are not valid, there is currently no warning in this case.
    '''
    auto_vectorize = False
    def __init__(self, wls = None, nis = None, itype='spline'):
        #Store RI if given, otherwise assume that wlr,wli,nir,nii exist
        if wls is not None:
            self.wlr = self.wli = array(wls)
            self.nir = real(nis)
            self.nii = imag(nis)

        #Wavelengths must be in acending order
        sorted_indicies = argsort(self.wlr)
        self.wlr = self.wlr[sorted_indicies]
        self.nir = self.nir[sorted_indicies]

        #Interpolate real & imaginary parts separately
        sorted_indicies = argsort(self.wli)
        self.wli = self.wli[sorted_indicies]
        self.nii = self.nii[sorted_indicies]

        self.itype = itype
        self.wl_limits = min(self.wlr.min(), self.wli.min()), max(self.wlr.max(), self.wli.max())
        
    def spline_interpolate(self, wls):
        from scipy import interpolate
        wls = atleast_1d(wls)
        
        ni = zeros(wls.shape, dtype=complex)
        for ii in range(len(ni)):
            #Constuct local spline interpolation Np points about wl
            Np = 4
            iselectr = absolute(self.wlr-wls[ii]).argmin()
            iselecti = absolute(self.wli-wls[ii]).argmin()
            wlr = self.wlr[max(iselectr-Np,0):iselectr+Np]
            nir = self.nir[max(iselectr-Np,0):iselectr+Np]
            wli = self.wli[max(iselecti-Np,0):iselecti+Np]
            nii = self.nii[max(iselecti-Np,0):iselecti+Np]
        
            #Evaluate local spline interpolation
            interpr = interpolate.UnivariateSpline(wlr, nir)
            interpi = interpolate.UnivariateSpline(wli, nii)
            
            ni[ii] = squeeze(interpr(wls[ii]) + interpi(wls[ii])*1j)
        return ni

    def linear_interpolate(self, wl):
        from scipy import interp
        xr = interp(atleast_1d(wl), self.wlr, self.nir)
        xi = interp(atleast_1d(wl), self.wli, self.nii)
        return squeeze(xr+1j*xi)

    def index_function(self, wavelength):
        if self.itype=='spline':
            return self.spline_interpolate(wavelength)
        else:
            return self.linear_interpolate(wavelength)

    def plot(self,wlrange=None, points=500, showdata=False):
        '''
        Plot the material refractive index
        * wrange is the wavelength range
        * points is the number of points to use in the plot
        '''
        if showdata:
            import Plotter as pl
            if wlrange:
                selectr = nonzero((min(wlrange)<self.wlr) & (self.wlr<max(wlrange)))
                selecti = nonzero((min(wlrange)<self.wli) & (self.wli<max(wlrange)))
            else:
                selectr = selecti = slice(None)
            
            pl.plot(self.wlr[selectr], self.nir[selectr],'x')
            pl.plot(self.wli[selecti], self.nii[selecti],'.')

        #Do default plot as well
        Material.plot(self,wlrange,points)
        

class GeneralFile(Interpolated):
    '''
    Refractive index from file formatted in columns of
    wavelength      nreal   nimag
    '''
    def __init__(self, mfile, wlunits='um'):
        #Check if we have a file or filename
        if hasattr(mfile,'readlines'):
            f = mfile
        #Otherwise try to open the file
        else:
            try:
                f = open(mfile, 'r')
            except:
                logging.error("Couldn't load material data file '%s'" % mfile)
                raise IOError
        
        #Read data from file
        data = self.readfile(f, 3)
        
        #Convert wls
        wls, nr, ni = data[:,:3].T
        wls = self.convertwavelength(wls, wlunits)
        
        Interpolated.__init__(self, wls = wls, nis = nr+ni*1j)
        
    def readfile(self, f, linelength=3, sep=' '):
        data = []
        linedata = []
        for l in f.readlines():
            #Try and read in a line of numbers
            try:
                #Replace the sep character with spaces to split on the whitespace
                linedata += map(float, l.replace(sep,' ').split())

                if len(linedata)>=linelength:
                    data += [linedata[:linelength+1]]
                    linedata = []
                
            #Skip this line if there is anything else
            except ValueError:
                continue
        
        return asarray(data)

class SopraFile(GeneralFile):
    '''
    Database n&k

    SOPRA provides with their instruments one of the largest database of
    optical indices available in the world (278 materials ).  We have now
    decided to give this database for free to the scientific community.
    It can be download very easily here.  Nevertheless the inverse can
    also be true and if you send us your own optical indices with the
    allowance to add them in the database their can be useful to other
    people.  This database comes with a help file which provide some more
    details about these data and also their origin.

    File Format:
    The different NK files are in ASCII format (*.NK files).  The first
    number indicates the spectral unit (1 = eV, 2 = <B5>m, 3 = cm-1, 4 = nm).
    The second and the third ones indicate the spectral range and the
    fourth one the number of intervals.  The N and K values are then
    stored in the order N,K,N,K... One example is reported below:

    1, 0.6, 5.9, 53
    3.4471, 0.
    3.4595, 0.
    ...
    1.1328, 3.0447
    1.0831, 2.9823
    '''
    def __init__(self, mfile):
        #Check if we have a file or filename
        if hasattr(mfile,'readlines'):
            f = mfile
        #Otherwise try to open the file
        else:
            try:
                f = open(mfile, 'r')
            except:
                logging.error("Couldn't load material data file '%s'" % mfile)
                raise IOError

        #Read in first info line
        try:
            params = map(float, f.readline().split())
            if len(params)<>4: raise ValueError
        except ValueError:
            raise ValueError, "File %s not a SOPRA data file!"
        
        wlunits = {1:'ev', 2:'um', 3:'/cm', 4:'nm'}[int(params[0])]
        
        #Reconstruct the wls
        dx = (params[2]-params[1])/params[3]
        datarange = arange(params[1],params[2]+dx,dx)
        
        #Convert wls to um
        wls = self.convertwavelength(datarange, wlunits)
        
        #Read data from file
        data = self.readfile(f, 2)
        nr = data[:,0]; ni = data[:,1]
        
        Interpolated.__init__(self, wls = wls, nis = nr+ni*1j)


# ------------------------- SOPRA interpolated data -----------------------------

class Sopra(SopraFile):
    def __init__(self, name):
        self.name = "Sopra data from %s" % name
        import pkg_resources
        
        respath = 'data/Sopra/'
        if pkg_resources.resource_exists(__name__, respath+'%s.all' % name):
            raise NotImplementedError, "Multi-file luxpop data not implemented"
        elif pkg_resources.resource_exists(__name__, respath+'%s.nk' % name):
            sopra_file = pkg_resources.resource_stream(__name__, respath+'%s.nk' % (name))
        else:
            raise RuntimeError, "Couldn't find Sopra datafile"
        
        SopraFile.__init__(self,sopra_file)

# ------------------------- LUXPOP interpolated data -----------------------------

class Luxpop(GeneralFile):
    """
    Load a material from the Luxpop database (not included)
    it can be downloaded from http://www.luxpop.com/RefractiveIndexList.html
    and put into the directory "Polymode/data/luxpop" before installing
    """
    def __init__(self, names):
        name = names
        self.name = "Luxpop data from %s" % name
        
        #Try and locate file
        from pkg_resources import resource_stream, resource_exists
        
        respath = 'data/luxpop/'
        
        #Try and find file with the different formats
        if resource_exists(__name__, respath+'%s.nk' % name):
            wlunits = "A"
            luxpop_file = resource_stream(__name__, respath+'%s.nk' % name)
            data = self.readfile(luxpop_file ,3)
            wlr = wli = data[:,0]; nreal = data[:,1]; nimag = data[:,2]

        elif resource_exists(__name__, respath+'%s_n.dat' % name):
            wlunits = "nm"
            luxpop_filen = resource_stream(__name__, respath+'%s_n.dat' % name)
            datan = self.readfile(luxpop_filen, 2, sep=',')
            wlr = datan[:,0]; nreal = datan[:,1]

            #Use _k file for imaginary part
            if resource_exists(__name__, respath+'%s_k.dat' % name):
                luxpop_filek = resource_stream(__name__, respath+'%s_k.dat' % name)             
                datak = self.readfile(luxpop_filek, 2, sep=',')
                wli = datak[:,0]; nimag = datak[:,1]
            else:
                wli = wlr; nimag = nreal*0

        else:
            raise ValueError, "Material file not found"

        #Convert the wavelengths and store the data
        self.wlr = self.convertwavelength(wlr, wlunits)
        self.wli = self.convertwavelength(wli, wlunits)
        self.nir, self.nii = nreal,nimag

        Interpolated.__init__(self)

# ------------------------- nk interpolated data from FreeSnell -----------------------------

class FreeSnell(GeneralFile):
    """
    Load a material from the internal Freesnell materials database
    This is taken from the freesnell thinfilm optical simulator at
    http://people.csail.mit.edu/jaffer/FreeSnell/
    """
    def __init__(self, filename=None):
        from pkg_resources import resource_stream
        filename = self.filename if filename is None else filename
        freesnell_file = resource_stream(__name__, 'data/freesnell/%s.nk' % filename)
        GeneralFile.__init__(self, freesnell_file, 'ev')

class Iron(FreeSnell):
    name = "Iron (Freesnell)"
    citation = "Weaver, J. H., Colavita, E., Lynch, D. W., and Rosei, R., Phys. Rev. Sect. B, 19, 3850, 1979."
    color = '#B7410E'
    filename = 'fe'

class Chromium(FreeSnell):
    name = "Chromium (Freesnell)"
    color = ''
    citation = "Bos, L. W., and Lynch, D. W., Phys. Rev. Sect. B, 2, 4567, 1970."
    filename = 'cr'

class Silver(FreeSnell):
    name = "Silver (Freesnell)"
    citation = "Hagemann, H. J., Gudat, W., and Kunz, C., J. Opt. Soc. Am., 65, 742, 1975."
    color = 'silver'
    filename = 'ag'

class Gold(FreeSnell):
    name = "Gold (Freesnell)"
    citation = "Olson, C. G., Lynch, D. W., and Weaver, J. H., unpublished."
    color = 'gold'
    filename = 'au'

class Aluminium(FreeSnell):
    name = "Aluminium (Freesnell)"
    citation = "Shiles, E., Sasaki, T., Inokuti, M., and Smith, D. Y., Phys. Rev. Sect. B, 22, 1612, 1980."
    color = 'gainsboro'
    filename = 'al'

class Copper(FreeSnell):
    name = "Copper (Freesnell)"
    citation = "Hagemann, H. J., Gudat, W., and Kunz, C., J. Opt. Soc. Am., 65, 742, 1975."
    color = '#B87333'
    filename = 'cu'

class Nickel(FreeSnell):
    name = "Nickel (Freesnell)"
    citation = "Lynch, D. W., Rosei, R., and Weaver, J. H., Solid State Commun., 9, 2195, 1971."
    color = 'aliceblue'
    filename = 'ni'

class Lithium(FreeSnell):
    name = "Lithium (Freesnell)"
    citation = "Lynch, D. W., and Hunter, W. R., in HOC-II, p.345."
    color = 'aliceblue'
    filename = 'li'

class Vanadium(FreeSnell):
    name = "Vanadium (Freesnell)"
    citation = "Olson, C. G., Lynch, D. W., and Weaver, J. H., unpublished."
    color = 'aliceblue'
    filename = 'v'

class Tungsten(FreeSnell):
    name = "Tungsten (Freesnell)"
    citation = "Weaver, J. H., Lynch, D. W., and Olson, C. G., Phys. Rev. Sect. B, 12, 1293, 1975."
    color = 'aliceblue'
    filename = 'w'

class Zirconium(FreeSnell):
    name = "Zirconium, polycrystaline (Freesnell)"
    citation = "Lanham, A. P., and Terherne, D. M., Proc. Phys. Soc., 83, 1059, 1964."
    color = 'aliceblue'
    filename = 'zr'

class HDPE(FreeSnell):
    name = "High Density Polyethylene (Freesnell)"
    color = 'peru'
    filename = 'hppe'

class Germanium(FreeSnell):
    name = "Germanium, single crystal (Freesnell)"
    citation = "Potter, R. F., in HOC-I, p.465."
    color = 'palegreen'
    filename = 'ge'

class Silicon(FreeSnell):
    name = "Silicon, single crystal (Freesnell)"
    citation = "Edwards, D. F., in HOC-I, p. 547."
    color = 'mistyrose'
    filename = 'si'

class Water(GeneralFile):
    name = "Water (Stegelstein)"
    citation = "Segelstein, D., 1981 M.S. Thesis, University of Missouri--Kansas City"
    color = 'aqua'
    def __init__(self):
        from pkg_resources import resource_stream
        freesnell_file = resource_stream(__name__, 'data/other/water_segelstein.nk')
        GeneralFile.__init__(self, freesnell_file, 'um')

class Water_IR(GeneralFile):
    name = "Water (Wieliczka)"
    citation = "Wieliczka, Weng & Querry (Appl. Opt. 28, 1714-1719, 1989)"
    color = 'aqua'
    def __init__(self):
        from pkg_resources import resource_stream
        freesnell_file = resource_stream(__name__, 'data/other/water_wieliczka.nk')
        GeneralFile.__init__(self, freesnell_file, '/cm')


#==================== Fixed Index ========================#

class Air(Fixed):
    __doc__ = name = "Air"
    color = 'azure'
    def __init__(self):
        Fixed.__init__(self, 1.0)

#==================== Optical Polymers ========================#

class PMMA(CauchyApproximation):
    __doc__ = utf8out(u"Polymethyl methacrylate (C₅O₂H₈)n")
    name = "PMMA"
    color = 'dodgerblue'
    wl_limits = (0.25,1.5)
    A=[2.18645820,-2.4475348e-04,1.4155787e-02,-4.4329781e-04,7.7664259e-05,-2.9936382e-06]

Polymer = PMMA

class PMMAKasarova(CauchyApproximation):
    __doc__ = utf8out(u"Polymethyl methacrylate (C₅O₂H₈)n Polymer at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "PMMA"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A = [2.399964,-8.308636e-2,-1.919569e-1,8.720608e-2,-1.666411e-2,1.169519E-3]

class Polystyrene(CauchyApproximation):
    __doc__ = utf8out(u"Polystyrene (C8H9)n at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Polystyrene"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.610025, -6.143673E-2, -1.312267E-1, 6.865432E-2, -1.295968E-2, 9.055861E-4]

class Polycarbonate(CauchyApproximation):
    __doc__ = utf8out(u"Polycarbonate at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Polycarbonate"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.633127, -7.937823E-2, -1.734506E-1, 8.609268E-2, -1.617892E-2, 1.128933E-3]

class SAN(CauchyApproximation):
    __doc__ = utf8out(u"Styrene acrylonitrile at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "SAN"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.595568, -6.848245E-2, -1.459074E-1, 7.329172E-2, -1.372433E-2, 9.426682E-4]

class CTERich(CauchyApproximation):
    __doc__ = utf8out(u"CTE Rich.® at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "CTE Rich"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.663794, -1.059116E-1, -2.492271E-1, 1.165541E-1, -2.211611E-2, 1.545711E-3]

class NAS21(CauchyApproximation):
    __doc__ = utf8out(u"Methyl methacrylate styrene copolymer at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "NAS21"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.054612, 1.374019E-1, 3.200690E-1, -1.152867E-1, 2.077225E-2, -1.383569E-3]

class LowStyrene(CauchyApproximation):
    __doc__ = utf8out(u"Low Styrene® at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Low Styrene"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.360004, -4.014429E-2, -8.371568E-2, 4.160019E-2, -7.586052E-3, 5.071533E-4]

class Optorez(CauchyApproximation):
    __doc__ = utf8out(u"Optorez 1330® at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Optorez"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.291142, -3.311944E-2, -1.630099E-2, 7.265983E-3, -6.806145E-4, 1.960732E-5]

class Zeonex(CauchyApproximation):
    __doc__ = utf8out(u"Zeonex E48R® at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Zeonex"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.482396, -6.959910E-2, -1.597726E-1, 7.383333E-2, -1.398485E-2, 9.728455E-4]

class Bayer(CauchyApproximation):
    __doc__ = utf8out(u"Bayer® polycarbonate at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Bayer"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.542676, -4.366727E-2, -8.196872E-2, 4.718432E-2, -8.892747E-3, 6.324010E-4]

class Cellulose(CauchyApproximation):
    __doc__ = utf8out(u"Cellulose (laboratory specimen) at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Cellulose"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.139790, -6.317682E-3, -5.920813E-3, 9.613514E-3, -1.967293E-3, 1.363793E-4]

class Polyacrylate(CauchyApproximation):
    __doc__ = utf8out(u"Polyacrylate (laboratory specimen) at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Polyacrylate"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.364830, -6.955268E-2, -1.356107E-1, 6.053900E-2, -1.166640E-2, 8.542615E-4]

class Styrene(CauchyApproximation):
    __doc__ = utf8out(u"Styrene (laboratory specimen) at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Styrene"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.274658, -5.700326E-3, -7.262838E-3, 1.233343E-2, -2.481307E-3, 1.784805E-4]

class PolycarbonateLS(CauchyApproximation):
    __doc__ = utf8out(u"Polycarbonate (laboratory specimen)  at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "PCLS"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.496875, -5.014035E-2, -4.188992E-2, 1.732175E-2, -1.240544E-3, -1.977750E-5]

class PolystyreneLS(CauchyApproximation):
    __doc__ = utf8out(u"Polystyrene (laboratory specimen) at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "PSLS"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[2.721609, -9.982812E-2, -2.518650E-1, 1.269202E-1, -2.549211E-2, 1.867696E-3]

class Acrylic(CauchyApproximation):
    __doc__ = utf8out(u"Acrylic (laboratory specimen) at 20°C")
    citation = "S. N. Kasarova et al., Optical Materials 29, 1481-1490, 2007"
    link = "doi:10.1016/j.optmat.2006.07.010"
    name = "Acrylic"
    color = 'dodgerblue'
    temperature = 20
    wl_limits = (0.4358,1.052)
    A=[1.866120, 2.085454E-1, 4.806770E-1, -1.840693E-1, 3.424849E-2, -2.340796E-3]


#================= Glasses ===================#

class Silica(Sellmeier):
    __doc__ = name = utf8out(u"Fused Silicon Oxide (SiO₂)")
    citation = "Malitson 1965, http://www.opticsinfobase.org/abstract.cfm?URI=josa-55-10-1205"
    color = 'yellow'
    wl_limits = (0.2,3.8)
    B = [0.6961663, 0.4079426, 0.8974794]
    L = [0.0684043, 0.1162414, 9.896161]

class Sapphire(Sellmeier):
    __doc__ = name = "Sapphire (Ordinary waves)"
    color = 'mediumblue'
    B = [1.43134930, 0.65054713, 5.3414021]
    L = sqrt([5.2799261e-3, 1.42382647e-2, 3.25017834e2])

class BK7(Sellmeier):
    __doc__ = name = "Borosilicate crown glass (BK7)"
    color = 'lightcoral'
    B = [1.03961212, 0.231792344, 1.01046945]
    L = sqrt([6.00069867e-3,2.00179144e-2,1.03560653e2])

class Quartz(Laurent):
    __doc__ = name = "Crystalline Quartz (Ordinary waves)"
    color = 'ivory'
    A = [-1.25900000e-02, 2.38490000, 1.07900000e-02, 1.65180000e-04, \
        -1.94741000e-06, 9.36476000e-08 ]

class Germania(ClassiusMorlotti):
    __doc__ = name = utf8out(u"Germanium Oxide (GeO₂)")
    citation = "Sunak, H.R.D.; Bastien, S.P., Photonics Technology Letters V1 N6 142-145, 1989"
    color = 'lawngreen'
    wl_limits = (0.6,1.8)
    
    B = aa_([0.2271125649, 0.1099158881, 0.1052670953])
    L = aa_([0.060928804, 0.1419148170, 10.86114943])

class SiO2GeO2(ClassiusMorlotti):
    __doc__ = utf8out(u"Silica (SiO₂) doped with a molar fraction f of Germania (GeO₂)")
    name = "Silica doped with Germania"
    citation = "Sunak, H.R.D.; Bastien, S.P., Photonics Technology Letters V1 N6 142-145, 1989"
    color = 'lightgreen'
    wl_limits = (0.6,1.8)
    
    B = aa_([0.2045154578, 0.06451676258, 0.1311583151])
    L = aa_([0.06130807320, 0.1108859848, 8.964441861])
    D = aa_([-0.1011783769, 0.1778934999, -0.164179581])
    def __init__(self, f=0):
        if f>0.2:
            logging.warning(utf8out(u"SiO₂/GeO₂: Dopant levels greater than 20% may be inaccurate") )
        self.f = f

class SiO2Fl(ClassiusMorlotti):
    __doc__ = utf8out(u"Silica (SiO₂) doped with a molar fraction f of Flourine")
    name = "Silica doped with Flourine"
    citation = "Sunak, H.R.D.; Bastien, S.P., Photonics Technology Letters V1 N6 142-145, 1989"
    color = 'cornsilk'
    wl_limits = (0.6,1.8)
    
    B = aa_([0.2045154578, 0.06451676258, 0.1311583151])
    L = aa_([0.06130807320, 0.1108859848, 8.964441861])
    D = aa_([-0.05413938039, -0.1788588824, -0.07445931332])
    def __init__(self, f=0):
        if f>0.02:
            logging.warning(utf8out(u"SiO₂Fl: Dopant levels greater than 2% may be inaccurate") )
        self.f = f


#------------------------------------- Finding Materials ---------------------------------------#

#Find available materials
__doc__ = "Material: Classes for modelling materials"
__doc__ += "\n\n**** Inbuilt Materials ****\n"
allinstances = filter( lambda item: hasattr(item[1],'name'), locals().items() )
__doc__ += "\n".join( [ "%s:  %s" % (name, mat.name) for name,mat in allinstances ] )
__doc__ += "\n\nUse the Material.find_materials() function to find other materials"

def find_materials(search_string):
    "Search the materials datafiles for a material containing the search term"
    try:
        import pkg_resources
    except ImportError:
        pkg_resources = None

    #Search in inbuilt materials
    found = {'Inbuilt': []}
    ss = search_string.lower()
    for name,mat in allinstances:
        #Strip out module name
        mname = str(mat)[str(mat).rfind('.')+1:]

        #Search in code name
        if name.lower().find(ss)>=0:
            found['Inbuilt'] += [ "%s:\t%s" % (mname, mat.name) ]
        #Search in long name
        elif mat.name.lower().find(ss)>=0:
            found['Inbuilt'] += [ "%s:\t%s" % (mname, mat.name) ]
    
    #Search in data directory
    if pkg_resources is not None:
        for datadir in material_resource_dirs:
            if pkg_resources.resource_isdir(__name__, 'data/'+datadir):
                names = pkg_resources.resource_listdir(__name__, 'data/'+datadir)
                names = filter(lambda name: name.lower().find(ss)>=0, names)
                names = filter(lambda name: name.lower()[-4:]<>'html', names)
                found[datadir.capitalize()] = names

    print "\nMaterials found in internal datafiles:"
    for resource in found:
        if len(found[resource])>0:
            print "\n%s materials:" % resource
            print "\n".join(found[resource]).replace('.nk','').replace('.dat','').replace('.all','')
    print
    
