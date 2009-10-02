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
'''
Base class for the storage and manipulation of mode data

'''

from __future__ import division
import sys, time, logging

from numpy import *
from numpy.lib.scimath import sqrt

from .mathlink  import utf8out, misc, timer, coordinates, constants, hankel1, hankel1p, jv
from .difflounge import finitedifference
from . import Plotter, Waveguide

def branchsqrt(x):
    '''
    Take the correct branch of the sqrt.
    The branchcut is along the -ve imaginary axis
    see: M. Nevière, "Electrodynamic theory of gratings", Chap 5., 1980.
    '''
    argx = remainder(angle(x)+pi/2,2*pi)-pi/2
    f = absolute(x)**0.5*exp(0.5j*argx)
    return f

def swap_vector(v, axis=-1, shift=0):
    "invert fourier frequencies"
    v = v.swapaxes(axis,0)
    if shift:
        #Flip the +/- frequency components in standard fft ordering
        fftflip = fft.fftshift(v,axes=[0])[::-1]

        #We need to roll the array by 1 is it's even
        if mod(v.shape[0],2)==1:
            v = fft.ifftshift(fftflip,axes=[0])
        else:
            v = fft.fftshift(roll(fftflip,1,axis=0),axes=[0])
    else:
        #Standard flip
        v = v[::-1]
    v = v.swapaxes(0,axis)
    return v

# +-----------------------------------------------------------------------+
# | Utility Functions
# +-----------------------------------------------------------------------+

def compress_modes(modes, Nshape=None, coord=None):
    """
    Compress modes to given shape or coord
    """
    for m in modes:
        if coord is None:
            coord = m.coord.new(Nshape=Nshape)
        m.compress(coord)
        m.add_extensions()
    return modes

def filter_suprious_modes(modes):
    """Select those modes that are converged and not spurious"""
    notspurious = nonzero([not md.is_spurious for md in modes])[0]
    converged = nonzero([md.is_converged() for md in modes])[0]
    wanted = intersect1d(notspurious, converged)
    
    modes = [modes[ii] for ii in wanted]
    return modes

def filter_unique_modes(modes, cutoff=1e-8):
    """Select those modes that are unique"""
    unique = ones(len(modes))
    smodes = sorted(modes)

    for ii in range(len(smodes)-1):
        m1 = smodes[ii]
        m2 = smodes[ii+1]
        if abs(m1.evalue-m2.evalue)<abs(m1.evalue)*cutoff:
            ip = dot(conj(m1.left), m2.right)
            if abs(ip)>1e-8:
                print "Repeated:",m1.neff,m2.neff
                unique[ii+1] = 0
    
    modes = [smodes[ii] for ii in transpose(nonzero(unique))]
    logging.info("Filtered %d repeated modes" % (len(smodes)-len(modes)))
    return modes

def construct_degenerate_pair(m1):
    """Constrcut degenerate mode pair from circularly polarized mode"""
    if m1.m0==0:
        logging.warning("Mode class 0 is non-degenerate.")

    #Construct opposite circularly polarized mode
    mode_class = type(m1)
    m2 = mode_class(coord=m1.coord, m0=-m1.m0, wl=m1.wl, evalue=m1.evalue)
    m2.exterior_index = m1.exterior_index
    m2.interior_index = m1.interior_index

    #Form new fields as hcalculated_electric_fielde+=conj(he-), he+=conj(he-)
    rs = conj(m1.shape_vector(m1.right)[...,::-1])
    ls = conj(m1.shape_vector(m1.left)[...,::-1])

    #Conjugation flips fourier coefficients
    m2.right = swap_vector(rs,1,1).ravel()
    m2.left = -swap_vector(ls,1,1).ravel()

    #Remember to add extensions
    m2.add_extensions()

    return m1,m2
    
def construct_lp_degenerate_pair(m):
    """Construct linearly polarized mode pair"""
    m1c,m2c = construct_degenerate_pair(m)

    #Construct LP modes as a combined vector
    # so m1 =(m1c+m2c)/2
    #     m2 =(m1c-m2c)/2
    m1 = CombinedVectorMode((m1c,m2c), (0.5,0.5), coord=m.coord, m0=0, 
                wl=m.wl, evalue=m.evalue)
    m2 = CombinedVectorMode((m1c,m2c), (-0.5j,0.5j), coord=m.coord, m0=0, 
                wl=m.wl, evalue=m.evalue)

    #remember to set these
    m1.exterior_index = m2.exterior_index = m.exterior_index
    m1.interior_index = m2.interior_index = m.interior_index
    
    return m1,m2

def construct_combined_mode(modes, coeffs, neff=None, wl=None):
    """
    Construct a mode from a linear combination of other modes.
    given the modes and the coefficients the mode fields are
    combined as:
        F = Σᵢ aᵢ Fᵢ
    Where F is the field Fᵢ is the field of the ith mode and aᵢ the
    ith coefficient.
    """
    assert len(modes)==len(coeffs), "Number of modes and coefficients must be equal"
    assert len(modes)>0, "At least one mode is required"
     
    #Take a representitive mode
    m = modes[0]
    coord = m.coord
    ev = (neff*m.k0)**2 if neff else m.evalue
    wl = wl if wl else m.wl
    
    mc = CombinedVectorMode(modes, coeffs, coord=coord, m0=0, wl=wl, evalue=ev)

    #remember to set these
    mc.exterior_index = m.exterior_index
    mc.interior_index = m.interior_index
    
    return mc


# +-----------------------------------------------------------------------+
# |
# | The Mode class
# |
# +-----------------------------------------------------------------------+
class Mode(object):
    """
    Base class for a mode, to construct a mode all paramters are optional
    and the mdoe information can be added later.
    
    Paramters:
     coord: internal coordinate object
     wl: the wavelength of the mode solution
     m0: the mode class
     evalue: the eigenvalue, automatically converted to beta and neff
     wg: the waveguide used to solve. not stored, just used to extract info
     left: the left eigenvector
     right: the right eigenvector
    """
    def __init__(self, coord=None, wl=1, m0=1, evalue=0, wg=None,
            symmetry=1, left=None, right=None):
        self.m0 = m0
        self.evalue = evalue
        self.wl=wl
        self.coord = coord
        self.symmetry = symmetry
        self.left = left
        self.right = right
        
        #Solver paramters
        self.store_field = True
        self.tolerance = 1e-9
        
        if wg:
            self.exterior_index = wg.exterior_index(wl)
            self.interior_index = wg.interior_index(wl)

        #Analytic extensions
        self.aleft_ext = self.aright_ext = None
        self.aleft_int = self.aright_int = None
        
        #Information
        self.convergence=[]
        self.is_spurious = False
        self.residue=0
        self.label = {}

    def copy(self, link_fields=False):
        """
        Create new copy of mode.
        link_vectors: share the mode information between modes
            to save memory
        """
        newmode = self.__class__()
        newmode.m0 = self.m0
        newmode.evalue = self.evalue
        newmode.wl = self.wl
        newmode.coord = self.coord
        newmode.symmetry = self.symmetry

        if hasattr(self,'exterior_index'):
            newmode.exterior_index = self.exterior_index
            newmode.interior_index = self.interior_index

        #Information
        newmode.tolerance = self.tolerance
        newmode.convergence = self.convergence
        newmode.is_spurious = self.is_spurious
        newmode.residue = self.residue
        newmode.label = self.label
        return newmode

    def floquet(self, m0=None, coord=None):
        m0 = self.m0 if m0 is None else m0

        if coord is None:
            return exp(1j*m0*self.coord.phiv)
        else:
            return exp(1j*m0*coord.phiv)

    # -------------------------------------------------------------------
    # Numerical properties of the Mode
    # -------------------------------------------------------------------

    @property
    def gamma_ext(self):
        """
        Calculate the mode parameter, ɣ = √[n₀² k₀² - β²]
        """
        n = self.exterior_index
        return branchsqrt((n*self.k0)**2-self.evalue)

    @property
    def gamma_int(self):
        """
        Calculate the mode parameter, ɣ = √[n₀² k₀² - β²]
        """
        n = self.interior_index
        return branchsqrt((n*self.k0)**2-self.evalue)

    @property
    def k0(self):
        "Return the wavenumber"
        return 2*pi/self.wl

    #Calulate/set neff from/to internal evalue
    # *** Note - setter isn't implemented in python 2.5 ***
    #@property
    def get_neff(self):
        "Return the effective index of the mode"
        return sqrt(self.evalue)/self.k0
    #@neff.setter
    def set_neff(self, value):
        "Return the effective index of the mode"
        self.evalue = value**2*self.k0**2
    neff = property(get_neff, set_neff)
    
    #@property
    def get_beta(self):
        "Return the modal eigenvalue, β"
        return sqrt(self.evalue)
    #@beta.setter
    def set_beta(self, value):
        "Return the modal eigenvalue, β"
        self.evalue =  value**2
    beta = property(get_beta, set_beta)
    
    @property
    def loss(self):
        "Calculate the confinement loss of the mode as loss = 2×10⁷ Im{β}/ln(10) db/m"
        return 2e7*imag(self.beta)/log(10)

    def is_converged(self):
        return self.residue<self.tolerance


    # -------------------------------------------------------------------
    # Mode information
    # -------------------------------------------------------------------
    def estimate_class(self, wg=None, core_size=None):
        "Guess the core mode class as either TM, TE or Hybrid"
        cutoff = 1e6
        cutoff_like = 100
        if core_size is None:
            core_size = wg.core_size*0.75
        #ec = self.coord.new( rmax=core_size, Nshape=(100,20), border=1 )
        ec = coordinates.PolarCoord(rrange=(0,core_size), arange=(-pi,pi), N=(50,50), border=1)
        
        #Look at the mode fields in the core
        h = self.magnetic_field(coord=ec)
        e = self.electric_field(wg, coord=ec)

        modeclass = 'HE/EH_(%d,j)-like' % mod(self.m0, self.symmetry)

        #All mode classes <>0 or symmetry/2 are degenerate
        if self.m0==0 or self.m0==self.symmetry/2.0:
            tm_factor = sqrt(abs(h[0])**2+abs(h[1])**2).sum()/median(abs(h[2]))
            te_factor = sqrt(abs(e[0])**2+abs(e[1])**2).sum()/median(abs(e[2]))
        
            if tm_factor/te_factor > cutoff:
                modeclass = 'TM_(0,j)'
            elif tm_factor/te_factor > cutoff_like:
                modeclass = 'TM_(0,j)-like'
            elif te_factor/tm_factor > cutoff:
                modeclass = 'TE_(0,j)'
            elif te_factor/tm_factor > cutoff_like:
                modeclass = 'TE_(0,j)-like'
        else:
            modeclass += " (degenerate)"
        return modeclass

    def polarization_angle(self, wg=None, core_size=None):
        "Approximate polarization angle of the mode (if linearly polarized)"
        if core_size is None:
            if wg is None:
                core_size = 0.75*self.coord.rmax
            else:
                core_size = wg.core_size*0.75

        ec = coordinates.PolarCoord(rrange=(0,core_size), arange=(-pi,pi), N=(50,50), border=1)
        
        #Look at the mode fields in the core
        hx,hy,hz = self.magnetic_field(cartesian=1, coord=ec)

        xpol = linalg.norm(hx)
        ypol = linalg.norm(hy)
        
        return arctan2(xpol,ypol)

    # -------------------------------------------------------------------
    # Overloadable methods
    # -------------------------------------------------------------------
    def discard_fields():
        pass
    def store_calculated_electric_field(self, wg=None, force=False):
        pass
    def add_extensions(self, bc_ext=-3, bc_int=2):
        pass
    def compress(self, size=None, astype=complex64):
        pass
    def normalize(self, wg=None, coord=None):
        pass

    # -------------------------------------------------------------------
    # Field properties
    # -------------------------------------------------------------------
    def poynting(self, **kwargs):
        '''
        The logitudinal component of the time averaged Poyting vector
        S_z = 0.5 ẑ∙(e×h*)
        '''
        ht = self.magnetic_transverse_field(**kwargs)
        et = self.electric_transverse_field(**kwargs)
        Sz = 0.5*cross(et,conj(ht),axis=0)
        return Sz

    def _convert_polar_vector(self, v, coord, cartesian=None):
        if cartesian is None: cartesian = not coord.polar_coordinate
        if cartesian: v = coordinates.vector_polar_to_cartesian(v, coord)       
        return v

    def mode_power(self, r=None, coord=None):
        '''
        The power in the computational region
        S_z = 1/2 |Aj|² ∫ (e×h*)∙ẑ dA
        '''
        if coord is None: coord = self.coord
        
        #Take a new maximum radius
        if r is not None and r<self.coord.rmax:
            Nshapenew = (int(r*coord.Nr/coord.rmax), coord.Naz)
            coord = coord.new(rmax=r, Nshape=Nshapenew)
        
        Sint = coord.int_dA(self.poynting(coord=coord))
        return Sint

    def mode_unconjugated_integral(self, fourier=True, r=None, extend=0, coord=None):
        '''
        I = ∫ (e×h)∙ẑ dA
        '''
        if coord is None: coord = self.coord
        
        #Take a new maximum radius
        if r is not None and r<self.coord.rmax:
            Nshapenew = (int(r*coord.Nr/coord.rmax), coord.Naz)
            coord = coord.new(rmax=r, Nshape=Nshapenew)
        
        ht = self.magnetic_transverse_field(fourier=fourier, coord=coord)
        et = self.electric_transverse_field(fourier=fourier, coord=coord)

        Sint = coord.int_dA(coord.cross_t(et,ht))
        
        #Add contibution from external fields
        #Does not accound for coord.symmetry currently - fix this!
        Sext = 0
        if extend or r>self.coord.rmax:
            Sext = -0.5j*iprod_extend(self.aleft_ext, self.aright_ext)
        if r>self.coord.rmax:
            Sext -= -0.5j*iprod_extend(self.aleft_ext, self.aright_ext, rc=r)
        
        if fourier: Sint/=coord.dphi    

        Sint = Sint*coord.symmetry + Sext*self.coord.dphi*self.symmetry
        return Sint

    def effective_area(self, coord=None):
        '''
        Nonlinear effective area of mode in um²
        Aeff = [∫ |Sz|² dA]² / ∫ |Sz|⁴ dA
        See Agrawal pg.
        '''
        if coord is None:   coord = self.coord
        
        Sz2 = real(self.poynting(coord=coord))**2
        Aeff = coord.int_dA(Sz2)**2 / coord.int_dA(Sz2**2)
        return Aeff
    
    def numerical_aperture(self, coord=None):
        '''
        Numerical aperture from the Gaussian approximation for a single moded MOF.
        See" Mortensen et al, "Numerical Aperture of a Single Mode Photonic Crystal Fiber"
        PTL, Vol. 14, No. 8, 2002, pp 1094-1096
        '''
        NA = 1/sqrt(1+pi*self.effective_area(coord=coord)/self.wl**2)
        return NA

    def spot_size(self, coord=None):
        '''
        Caculate the Petersen II spot size calculated as
        spot size = [∫ Sz² dA]² / ∫ (∇Sz)² dA
        '''
        if coord is None: coord=self.coord
        
        Sz = self.poynting(coord=coord)
        DSzr, DSzphi = coord.grad_t(Sz, m0=0, fourier=False) 
        p2 = coord.int_dA(Sz**2)/coord.int_dA(DSzr**2 + DSzphi**2)
        return sqrt(p2)

    def phase_velocity(self):
        '''
        Caculate the phase velocity (m/s) of the mode
        '''
        return phsycon.c*self.k0/self.beta

    def group_velocity(self, wg=None, coord=None):
        return constants.c/self.group_index(wg,coord)

    def group_index(self, wg=None, coord=None):
        '''
        Caculate the group index (m/s) of the mode based on the profile
        See Snyder and Love pg. 608
        '''
        if wg is None:      
                            logging.warning("Calculation being approximated as no waveguide available")
                            n2 = self.interior_index**2
                            dn2dl = 0
        else:
            n2 = wg.index2(wl=self.wl, coord=self.coord, resample=coord)
            dn2dl = wg.material.wavelength_derivative(self.wl, units='wl')

        h = self.magnetic_field(fourier=False, coord=coord)
        e = self.electric_field(wg, fourier=False, coord=coord)
        exh = cross(e[:2], conj(h[:2]), axis=0)
        e2abs = e*conj(e); h2abs = h*conj(h)
        
        if coord is None: coord=self.coord
        ng = coord.int_dA(h2abs + (n2-self.wl*dn2dl)*e2abs)/coord.int_dA(exh)/2
        return ng

    def integral_propagation(self, wg=None, coord=None):
        '''
        Caculate the propagation constant beta of the mode based on the
        mode profile, requires the waveguide
        '''
        if wg is None:      
                            logging.warning("Calculation being approximated as no waveguide available")
                            n2 = self.interior_index**2
        else:
            n2 = wg.index2(wl=self.wl, coord=self.coord, resample=coord)

        hr,hphi,hz = self.magnetic_field(fourier=0, coord=coord)
        er,ephi,ez = self.electric_field(wg, fourier=0, coord=coord)
        exh = cross((er,ephi), conj((hr,hphi)), axis=0)

        e2abs = er*conj(er)+ephi*conj(ephi)+ez*conj(ez)
        h2abs = hr*conj(hr)+hphi*conj(hphi)+hz*conj(hz)

        if coord is None: coord=self.coord
        beta = (2*self.k0)*coord.int_dA(n2*exh)/coord.int_dA(h2abs + conj(n2)*e2abs)
        return beta
        
    def integral_propagation_lossless(self, wg=None, coord=None):
        '''
        Caculate the propagation constant beta of the mode based on the
        mode profile, requires the waveguide
        '''
        if wg is None:      
                            logging.warning("Calculation being approximated as no waveguide available")
                            n2 = self.interior_index**2
        else:
            n2 = wg.index2(wl=self.wl, coord=self.coord, resample=coord)

        hr,hphi,hz = self.magnetic_field(coord=coord)
        er,ephi,ez = self.electric_field(wg, coord=coord)

        e2abs = er*conj(er)+ephi*conj(ephi)+ez*conj(ez)
        h2abs = hr*conj(hr)+hphi*conj(hphi)+hz*conj(hz)
        exh = er*conj(hphi)-ephi*conj(hr)


    # -------------------------------------------------------------------
    # Plotting
    # -------------------------------------------------------------------
    
    def plot(self, plottype='Sz', style='pcolor', part='real', cartesian=None, wg=None,
                Nx=None, sectors=None , rmin=None, rmax=None, cmap=None, coord=None, style1d='-',
                title=r"%(type)s, $n_{\mathrm{eff}}=%(tneff)s$"):
        """
        Plot the mode.
        
        Paramters:
         plottype: 'vector', 'Sz'*, 'Hz', 'Ez, 'Hr, 'Ha', 'Er, 'Ea', 'Ex', 'Ey', 'Hx', 'Hy'
         style: 'contour', 'pcolor'*, 'line'
         style1d: line style for 1d plotting
         part: one of 'real'*, 'imag', 'abs', 'phase' and optionally 'log'
         rmax: maxiumum plot radius
         Nx: plot sampling
         cmap: color map (from pylab.cm)
         wg: waveguide (for plotting Ez)
         title: custom title for the plot
        """
        plottype = plottype.lower()
        plotstyle={}
        color = None

        #if not self.has_mode_data():
        #   raise RuntimeError, "No vectors stored in this mode, cannot plot"

        symm = self.symmetry

        #If a coord is not specified we need to guess a good one
        if coord is None:
            #Resample output to a suitable size
            if Nx is None:
                if plottype.startswith('vect') or plottype in ['ph','pe']:
                    cartesian = 1
                    Nx = (20,20)
                else:
                    if cartesian:
                        Nx = (100,100)  #Default cartesian resolution
                    elif self.coord.shape is None or self.coord.shape[1]==1:
                        Nx = (200, 1)     #Keep radial coordintates
                    else:
                        Nx = (50, 100)   #Default polar resolution

            if rmax is None: rmax = self.coord.rrange[1]
            if rmin is None: rmin = 0
            dT = pi if sectors is None else sectors*pi/symm
        
            #New plotcoord for resampled plotting
            if cartesian:
                plotcoord = coordinates.CartesianCoord(X=rmax, Y=rmax, N=Nx)
            else:
                plotcoord = coordinates.PolarCoord(rrange=(rmin,rmax), arange=(-dT,dT), N=Nx, border=0)
        else:
            plotcoord = coord
            
        #Choose what to plot
        if plottype == 'sz' or plottype == 'power':
            plotdata = (self.poynting(coord=plotcoord))
        elif plottype == 'szc':
            plotdata = (self._calculated_poynting(wg=wg, coord=plotcoord))
        elif plottype == 'vector' or plottype == 'vectorh' or style == 'vector':
            plotdata = self.magnetic_transverse_field(cartesian=1, coord=plotcoord)
            color='black'; plottype = 'H-Vector'; style = 'vector'
        elif plottype == 'vectore':
            plotdata = self.electric_transverse_field(cartesian=1, coord=plotcoord)
            color='blue'; plottype = 'E-Vector'; style = 'vector'
        elif plottype == 'ph':
            plotdata = self.magnetic_transverse_field(cartesian=1, coord=plotcoord)
            color='black'; plottype = 'H-Polarization'; style = 'circ'; part='all'
        elif plottype == 'pe':
            plotdata = self.electric_transverse_field(cartesian=1, coord=plotcoord)
            plotstyle['ec']='blue'; plottype = 'E-Polarization'; style = 'circ'; part='all'
        elif plottype == 'hz':
            plotdata = self.magnetic_field(coord=plotcoord)[-1]
        elif plottype == 'ez':
            plotdata = self.electric_field(wg, coord=plotcoord)[-1]
        elif plottype == 'hr':
            plotdata = self.magnetic_transverse_field(cartesian=0, coord=plotcoord)[0]
        elif plottype == 'ha':
            plotdata = self.magnetic_transverse_field(cartesian=0, coord=plotcoord)[1]
        elif plottype == 'er':
            plotdata = self.electric_transverse_field(cartesian=0, coord=plotcoord)[0]
        elif plottype == 'ea':
            plotdata = self.electric_transverse_field(cartesian=0, coord=plotcoord)[1]
        elif plottype == 'erc':
            plotdata = self.calculated_electric_field(wg, coord=plotcoord)[0]
        elif plottype == 'eac':
            plotdata = self.calculated_electric_field(wg, coord=plotcoord)[1]
        elif plottype=='hx':
            plotdata = self.magnetic_transverse_field(cartesian=1, coord=plotcoord)[0]
        elif plottype=='hy':
            plotdata = self.magnetic_transverse_field(cartesian=1, coord=plotcoord)[1]
        elif plottype == 'ex':
            plotdata = self.electric_transverse_field(cartesian=1, coord=plotcoord)[0]
        elif plottype == 'ey':
            plotdata = self.electric_transverse_field(cartesian=1, coord=plotcoord)[1]

        else:
            logging.error("Plottype isn't recognised")
            return

        #Select real, imag or abs parts
        parts = part.replace(',',' ').split()
        if 'abs' in parts: plotdata = abs(plotdata)
        elif 'phase' in parts: plotdata = arctan2(real(plotdata),imag(plotdata))
        elif 'imag' in parts: plotdata = imag(plotdata)
        elif 'real' in parts: plotdata = real(plotdata)
        if 'log' in parts: plotdata = log(abs(plotdata))
        
        #2D or 1D plot
        if 'vector' in style or plotdata.shape[1]>1:
            Plotter.plot_v( plotcoord, plotdata, style=style, cmap=cmap, color=color)

        else:
            Plotter.plot(plotcoord.rv, plotdata.squeeze(), style1d)
            Plotter.xlabel('Radial distance, r')

        tdata = {'type':plottype, 'rneff':real(self.neff), 'ineff':imag(self.neff), 'loss':self.loss, 'm0':self.m0}
        tdata['tneff'] = misc.format_complex_latex(self.neff)
    
        Plotter.title( title % tdata )
        return plotcoord

    def info(self, wg=None):
        info_str = self.__str__().decode('utf8') + "\n"
        info_str += u" | Effective area: %.5g μm²\n" % self.effective_area()
        info_str += u" | Spot size: %.5g μm\n" %(self.spot_size())
        info_str += u" | Single moded numerical aperture: %.5g\n" %(self.numerical_aperture())
        if self.is_spurious:
            info_str += " | Possible spurious mode\n"

        if wg:
            quiet = wg.quiet; wg.quiet = True
            rc = wg.core_size

            info_str += " | Group index: %s m/s\n" % self.group_index(wg)
            info_str += " | Mode class: %s\n" % self.estimate_class(wg)
            info_str += " | Power in core: %.5g%%\n" % (100*self.mode_power(r=rc)/self.mode_power())

            wg.quiet = quiet
        print utf8out(info_str)

    # -------------------------------------------------------------------
    # Misc functions for 2 Modes
    # -------------------------------------------------------------------
    
    def __cmp__(self,other):
        "Compare two Modes, based on eigenvalue"
        if hasattr(other, "neff"):
            return cmp(self.neff,other.neff)
        else:
            return cmp(self.neff,other)

    def __mul__(self, other):
        "Construct numerical innerproduct as this.L dot that.R:"
        if hasattr(self,'right') and hasattr(other,'left'):
            return sum(self.right*other.left)
        else:
            raise IndexError, "Modes do not have left & right members!"

    def __str__(self):
        info_dict = {}

        info_dict['res'] = max(atleast_1d(self.residue))
        info_dict['userlab'] = ", ".join(["%s: %s" % (x,self.label[x]) for x in self.label])
        info_dict['m0'] = self.m0
        info_dict['wl'] = self.wl

        if self.coord is None:
            info_dict['shape'] = info_dict['symm'] = info_dict['rmax'] = "?"
        else:
            info_dict['shape'] = self.coord.shape
            info_dict['rmin'], info_dict['rmax'] = self.coord.polar_bounds()[:2]
            info_dict['symm'] = "C%d" % self.symmetry
        
        #Construct information string
        info_str = u"Mode, size: %(shape)s, symmetry: %(symm)s, m₀: %(m0)d\n" % info_dict
        info_str += u"λ: %(wl).4g, r: %(rmin).3g -> %(rmax).3g, res: %(res).2g\n" % info_dict
        info_str += u"neff=%s, loss=%.4gdB/m, %s" % (misc.format_complex(self.neff), self.loss, info_dict['userlab'])
        return utf8out(info_str)

    def __repr__(self):
        res = max(atleast_1d(self.residue))
        slab = ("","S")[self.is_spurious]
        clab = ("E","")[self.is_converged()]
        userlab = ", ".join(["%s: %s" % (x,self.label[x]) for x in self.label])

        info_str = u"<%s: m₀=%d λ=%.4g neff=%s r:%.2g %s [%s%s]>" \
            % (self.__class__.__name__, self.m0, self.wl, \
                misc.format_complex(self.neff), res, userlab, slab, clab)
        return utf8out(info_str)


# +-----------------------------------------------------------------------+
# | Mode class for the Scalar Wave Equation
# +-----------------------------------------------------------------------+

class ScalarMode(Mode):
    pass
    
# +-----------------------------------------------------------------------+
# | Mode class for the Fourier Decomposition Method
# +-----------------------------------------------------------------------+

#class FourierMode(Mode):
class VectorMode(Mode):
    pmax = 2
    r_axis=-2
    az_axis=-1
    pv = array([1, -1])

    reverse_left = False
    reverse_right = False
    reverse_p = False
    
    def copy(self, link_fields=False):
        """
        Create new copy of mode.
        link_vectors: share the mode information between modes
            to save memory
        """
        newmode = Mode.copy(self, link_fields)

        #Copy vectors or link them
        if link_fields or self.right is None:
            newmode.right = self.right
        else:
            newmode.right = self.right.copy()
        if link_fields or self.left is None:
            newmode.left = self.left
        else:
            newmode.left = self.left.copy()
        
        return newmode

    # +-----------------------------------------------------------------------+
    # | Pickling marshalling functions
    # | Compress the mode before saving if required
    # +-----------------------------------------------------------------------+

    def __getstate__(self):
        "Pickle all needed data, ignore cached data"
        state = self.__dict__.copy()
        ignore_list = ["Mdtn"]

        for ignore in ignore_list:
            if ignore in state:
                del state[ignore]

        if not self.store_field:
            state['left'] = state['right'] = None
        
        return state
    
    def __setstate__(self,state):
        "Restore pickled data"
        self.__dict__.update(state)

    # -------------------------------------------------------------------
    # Vector access routines
    # -------------------------------------------------------------------
    def add_extensions(self, bc_ext=-3, bc_int=2):
        """
        Set up the extensions to model the fields
         at locations below r_min and above r_max
         """
        if hasattr(self,'right') and self.right is not None:
            self.aright_ext = AnalyticExtension(bc_ext)
            self.aright_int = AnalyticExtension(bc_int, interior=1)

            self.aright_ext.set_extension_mode(self)
            self.aright_int.set_extension_mode(self)

        if hasattr(self,'left') and self.left is not None:
            self.aleft_ext = AnalyticExtensionLeft(bc_ext)
            self.aleft_int = AnalyticExtensionLeft(bc_int, interior=1)

            self.aleft_ext.set_extension_mode(self)
            self.aleft_int.set_extension_mode(self)

    def compress(self, wg, size, astype=complex64):
        "Replace the calculated vectors with resampled versions"
        if size is None or size==self.coord.shape:
            logging.warning("Compressing to same size as mode, doing nothing")
            return
        
        newcoord = wg.get_coord(size)
        newshape = (newcoord.Nr, newcoord.Naz, self.pmax)
        
        #Resample the right mode vector if there is one
        right = self.right
        left = self.left
        if right is not None:
            right = self.shape_vector(right).transpose((2,0,1))
            right = newcoord.fourier_resample(right, self.coord, fourier=True)
            right = self.unshape_vector(right.transpose((1,2,0))).astype(astype)

        #Resample the left mode vector if there is one
        if left is not None:
            left = self.shape_vector(left).transpose((2,0,1))
            left = newcoord.fourier_resample(left, self.coord, fourier=True)
            left = self.unshape_vector(left.transpose((1,2,0))).astype(astype)

        self.coord = newcoord
        self.left = left
        self.right = right          
        self.add_extensions()

    def get_right(self, swap=False):
        """
        Return the right eigenvector of the mode.
        swap: swap the Fourier frequencies m⟶-m
        """
        if self.right is None:
            raise RuntimeError, "No right eigenvector is available!"

        #Shape the vector correctly
        mdata = rollaxis(self.shape_vector(self.right),2)

        if swap:
            mdata = swap_vector(mdata, self.az_axis, 1)
        return mdata

    def get_left(self, field=True, swap=False):
        """
        Return the left eigenvector of the mode.
        field:   return mode corrected to electric field of mode
        swap: swap the Fourier frequencies m⟶-m
        """
        if self.left is None:
            raise RuntimeError, "No left eigenvector is available!"

        #Shape the vector correctly
        mdata = rollaxis(self.shape_vector(self.left),2)
        
        #'Fix' up left field as electric field
        # so E+ = L+/r, E- = -L-/r
        if field:
            mdata = mdata/self.coord.rv[:,newaxis]
            mdata *= self.pv[:,newaxis,newaxis]

        if self.reverse_p:
            mdata = conj(mdata)
            mdata = mdata[::-1]
        
        if swap:
            mdata = swap_vector(mdata, self.az_axis, 1)
        return mdata

    def transform_left(self, factor=1):
        if self.reverse_p:
            self.left *= conj(factor)
        else:
            self.left *= factor

    def transform_right(self, factor=1):
        self.right *= factor

    def shape_vector(self,v):
        "Return vector shaped as (r,phi,p)"
        shape = self.coord.shape+(2,)
        v = v.reshape(shape)
        return v
    def unshape_vector(self,v):
        "Return vector shaped for storage"
        return v.ravel()
    @property
    def shape(self):
        shape = self.coord.shape+(self.pmax,)
        return shape

    def discard_fields(self):
        "Discard all field information"
        self.left = None
        self.right = None
        self.aleft_int = self.aleft_ext = None
        self.aright_int = self.aright_ext = None

    def normalize(self, by='field', arg=0, wg=None, coord=None):
        """
        Normalize the fields so the electric and magnetic vectors have the correct
        correspondance, including the intrinsic impedence of free space:
        H = Htrue, E = (e0/mu0)^1/2 Etrue
        """
        if self.coord is None or self.right is None: return

        #Normalize largest component of magnetic transverse vector to be real
        #ht = self.magnetic_transverse_field(fourier=1)
        #self.right /= exp(1j*arg)*misc.absmax(ht)/abs(misc.absmax(ht))

        if self.left is None: return

        #Normalize E/H fields:
        P0 = self.mode_power(coord=coord)
        enorm = None
        if by=='poynting':
            enorm = conj(misc.absmax(self.poynting(coord=coord)))
        elif by=='ext':
            self.add_extensions()
            enorm = self._normalize_electric_field_extension()
        elif by=='field':
            enorm = self._normalize_electric_field(wg=wg,coord=coord)

        #Otherwise just use the power
        if enorm is None or isnan(enorm):
            enorm = 1./P0
        
        #Normalize absolute power so |P|=1
        #Can't normalize power phase and field relationship simultaneously
        b = sqrt(abs(P0*enorm))
        self.transform_right(1/b)
        self.transform_left(enorm/b)

        #Update analytic extensions
        self.add_extensions()

        #Recalulate power & field for information
        Pang = angle(self.mode_power(coord=coord))
        logging.debug(u"Normalized mode to power angle ∠%.3gπ" % (Pang/pi))
        return enorm
        
    def _normalize_electric_field(self, wg, fourier=True, coord=None):
        """
        Normalize Electric/Magnetic components so they have the correct relationship
        using the analytic extension
        H = Htrue, E = (e0/mu0)^1/2 Etrue
        """
        #Select the internal nodes only (electric vector will be incorrect at the boundary
        #nodes
        E = self.electric_transverse_field(fourier=fourier, coord=coord)[:,2:-2]
        Ec = self.calculated_electric_field(wg=wg, fourier=fourier, coord=coord)[:2,2:-2]

        xselect = absolute(E).argmax()
        enorm = array(Ec).flat[xselect]/array(E).flat[xselect]
        return enorm

    def _normalize_electric_field_extension(self):
        """
        Normalize Electric/Magnetic components so they have the correct relationship
        using the analytic extension
        H = Htrue, E = (e0/mu0)^1/2 Etrue
        """
        ap,am = self.aleft_ext.alphac
        bp,bm = self.aright_ext.alphac
        
        #The field coefficients should be such that this is one
        chi = (bp-bm)/(am+ap)/(1j*self.beta/self.k0)
        w = abs(bp-bm)**2+abs(am+ap)**2

        #Check that we aren't comparing two zero numbers
        if mean(w)<1e-10:
            print "Extension norm failed!"
            return None
            
        #Calculate the weighted average if more than one Fourier component
        if shape(w)[0]>1:
            #Select only low frequency components
            msl = absolute(self.aright_ext.msx).max(axis=0)<20
            bnorm = trapz(chi[msl]*w[msl])/trapz(w[msl])
        else:
            bnorm = chi
        return bnorm

    def _resample_vector(self, v, ext=(None,None), fourier=False, cartesian=None, coord=None):
        "Resample raw vectors to new coordinate"
        #Resample if coord is given and is not self.coord
        if coord and coord!=self.coord:
            v = coord.fourier_resample(v, self.coord, ext=ext, m0=self.m0, fourier=fourier)
        elif not fourier:
            v = fft.ifft(v, axis=self.az_axis)*self.floquet()
        
        #Convert to cartesian if specified
        if coord is None: coord=self.coord
        if cartesian is None: cartesian = not coord.polar_coordinate
        if cartesian: v = coordinates.vector_polar_to_cartesian(v, coord)       
        return v

    def calculated_electric_field(self, wg=None, fourier=False, cartesian=False, coord=None):
        "The transverse electric vector calculated from the internal H field"
        if wg is None:      
                            logging.warning("CEF Electric field calculation being approximated as no waveguide available")
                            n2 = self.interior_index**2
        else:
            n2 = wg.index2(wl=self.wl, coord=self.coord)

        hr,hphi,hz = self.magnetic_field(fourier=1)
        Dhzr,Dhzphi = self.coord.grad_t(hz, m0=self.m0, fourier=1)

        er = fft.ifft(self.beta*hphi+1j*Dhzphi, axis=1)/(self.k0*n2)
        ephi = -fft.ifft(self.beta*hr+1j*Dhzr, axis=1)/(self.k0*n2)
        ez = fft.ifft(1j*self.coord.curl_t((hr,hphi), fourier=1, m0=self.m0))/(self.k0*n2)
        e = array([er,ephi,ez])
        
        e = self._resample_vector(fft.fft(e,axis=-1), fourier=fourier, cartesian=cartesian, coord=coord)
        return e

    def store_calculated_electric_field(self, wg=None, force=False):
        '''
        Calculate the electric vector from the internal H field
        and save it as the internal electric vector
        '''
        if not (self.left is None or force):
            return
        logging.debug("Calculating electric field")
        
        if wg is None:      
            logging.warning("Electric field calculation being approximated as no waveguide available")
            n2 = self.interior_index**2
        else:
            n2 = wg.index2(wl=self.wl, coord=self.coord)
        
        hr,hphi,hz = self.magnetic_field(fourier=1)
        Dhzr,Dhzphi = self.coord.grad_t(hz, m0=self.m0, fourier=1)

        #Caclulate transverse electric field from Maxwell's equations
        er = fft.ifft(self.beta*hphi+1j*Dhzphi, axis=1)/(self.k0*n2)
        ephi = -fft.ifft(self.beta*hr+1j*Dhzr, axis=1)/(self.k0*n2)

        #Calculate vector
        ev = array([er+ephi*1j,er-ephi*1j])*self.coord.rv[:,newaxis]/self.pv[:,newaxis,newaxis]
        ev = fft.fft(ev,axis=-1)

        #Reverse if called for
        if self.reverse_left:
            ev = swap_vector(ev, -1, 1)
            
        self.left = self.unshape_vector(ev.transpose((1,2,0)))

    def magnetic_transverse_field(self, fourier=False, cartesian=None, coord=None):
        '''
        The transverse magnetic field, calculated from the internal H⁺,H⁻"
        cartesian=False returns h_t=(h_r,h_ϕ)
        cartesian=True returns h_t=(h_x,h_y)
        '''
        hp,hm = self.get_right(swap=self.reverse_right)
        ht = asarray([(hp+hm)/2, (hp-hm)/2j])

        ext = (self.aright_int, self.aright_ext)
        ht = self._resample_vector(ht, ext, fourier, cartesian, coord)
        return ht

    def magnetic_field(self, fourier=False, cartesian=None, coord=None):
        """
        The three component magnetic field (h_r,h_ϕ,h_z) or (hx,hy,hz)
        """
        hr,ha = self.magnetic_transverse_field(fourier=1)
        hz = 1j/self.beta * self.coord.div_t((hr,ha), fourier=1, m0=self.m0)
        h =asarray([hr,ha,hz])
        
        ext = (self.aright_int, self.aright_ext)
        h = self._resample_vector(h, ext, fourier, cartesian, coord)
        return h

    def electric_transverse_field(self, fourier=False, cartesian=None, coord=None, wg=None):
        '''
        The transverse electric field, calculated from the internal E⁺,E⁻"
        cartesian=False returns e_t=(e_r,e_ϕ)
        cartesian=True returns e_t=(e_x,e_y)
        
        if calculated_electric_field is true then calculate from the the magnetic field
        '''
        if self.left is None:
            #return self.calculated_electric_field(wg=wg, fourier=fourier, cartesian=cartesian, coord=coord)[:2]
            self.left = self._calculate_electric_vector(wg)

        ep,em = self.get_left(field=1, swap=self.reverse_left)
        et = asarray([(ep+em)/2, (ep-em)/2j])
        
        ext = (self.aleft_int, self.aleft_ext)
        et = self._resample_vector(et, ext, fourier, cartesian, coord)
        return et

    def electric_field(self, wg=None, fourier=False, cartesian=None, coord=None):
        """
        The three component electric field (e_r,e_ϕ,e_z)
        if calculated_electric_field is true then calculate from the the magnetic field
        """
        if self.left is None:
            #return self.calculated_electric_field(wg=wg, fourier=fourier, cartesian=cartesian, coord=coord)[:2]
            self.left = self._calculate_electric_vector(wg)

        if wg is None:
            logging.warning("Electric field calculation being approximated as no waveguide available")
            n2 = self.interior_index**2
        else:
            n2 = wg.index2(wl=self.wl, coord=self.coord)

        hr,ha = self.magnetic_transverse_field(fourier=1)
        er,ea = self.electric_transverse_field(fourier=1)
        ez = fft.fft(1j/(self.k0*n2) * fft.ifft(self.coord.curl_t((hr,ha), fourier=1, m0=self.m0), axis=1),axis=1)
        e = asarray([er,ea,ez])
        
        ext = (self.aleft_int, self.aleft_ext)
        e = self._resample_vector(e, ext, fourier, cartesian, coord)
        return e

# +-----------------------------------------------------------------------+
# 
# A cheap way to create linear combinations of modal fields
# without worrying about the technical details
# 
# +-----------------------------------------------------------------------+


class CombinedVectorMode(VectorMode):
    def __init__(self, modelist, coeffs, **kwargs):
        VectorMode.__init__(self, **kwargs)

        self.modes = modelist
        self.fieldcoeffs = coeffs
        
    def magnetic_transverse_field(self, *args, **kwargs):
        '''
        The transverse magnetic field, calculated from combined modes
        cartesian=False returns h_t=(h_r,h_ϕ)
        cartesian=True returns h_t=(h_x,h_y)
        '''
        ht=None
        for m,c in zip(self.modes, self.fieldcoeffs):
            if ht is None:
                ht = c*m.magnetic_transverse_field(*args, **kwargs)
            else:
                ht += c*m.magnetic_transverse_field(*args, **kwargs)
        return ht

    def electric_transverse_field(self, *args, **kwargs):
        '''
        The transverse electric field, calculated from the internal E⁺,E⁻"
        cartesian=False returns e_t=(e_r,e_ϕ)
        cartesian=True returns e_t=(e_x,e_y)
        '''
        et=None
        for m,c in zip(self.modes, self.fieldcoeffs):
            if et is None:
                et = c*m.electric_transverse_field(*args, **kwargs)
            else:
                et += c*m.electric_transverse_field(*args, **kwargs)
        return et

    def magnetic_field(self, *args, **kwargs):
        "The three component magnetic field (h_r,h_ϕ,h_z) or (hx,hy,hz)"
        h=None
        for m,c in zip(self.modes, self.fieldcoeffs):
            if h is None:
                h = c*m.magnetic_field(*args, **kwargs)
            else:
                h += c*m.magnetic_field(*args, **kwargs)
        return h

    def electric_field(self, *args, **kwargs):
        "The three component electric field (e_r,e_ϕ,e_z)"
        e=None
        for m,c in zip(self.modes, self.fieldcoeffs):
            if e is None:
                e = c*m.electric_field(*args, **kwargs)
            else:
                e += c*m.electric_field(*args, **kwargs)
        return e


# +-----------------------------------------------------------------------+
# 
#  Class for calulating the mode in the region beyond
#  the computational boundary.
# 
# +-----------------------------------------------------------------------+

class AnalyticExtension:
    '''
    Class for modeling the mode behaviour in the region beyond the computational boundary.
    '''
    def __init__(self, bcnode=-1, factor=1, interior=False):
        self.gamma = 0
        self.bc_node = bcnode
        self.interior = interior
        self.factor=factor

        if interior:
            self.bcoeffs = (0,0,1)
        else:
            self.bcoeffs = (1,0,0)

    def bessel(self, m, x):
        y = 0
        if self.bcoeffs[0]!=0:
            y += self.bcoeffs[0]*hankel1(m,x)
        if self.bcoeffs[1]!=0:
            y += self.bcoeffs[1]*hankel1(m,x)
        if self.bcoeffs[2]!=0:
            y += self.bcoeffs[2]*jv(m,x)
        
        #Make sure we don't get NaN's for large orders
        place(y, isnan(y), 0)
        return y
        
    def set_extension_mode(self, mode):
        self.factor = mode.beta
    
        if self.interior: gamma = mode.gamma_int
        else: gamma = mode.gamma_ext

        mdata = mode.get_right(swap=mode.reverse_right)
        return self.set_extension(mdata, mode.m0, mode.pv, mode.coord, gamma)

    def set_extension(self, mdata, m0, p, coord, gamma):
        self.rc = coord.rv[self.bc_node]
        self.ms, self.phi = coord.ms, coord.phiv
        self.m0, self.p = m0, p
        self.msx = self.ms + m0 + p[:,newaxis]

        #Take the points (in phi and p) at the correct r location
        #Note, need to copy it here, not have it linked to the main array
        self.mdata = mdata[...,self.bc_node,:] + 0
        self.gamma = gamma
        return self
        
    #Calculate analytic coefficients using the value of the mode at the boundary
    def calc_alphac(self, gamma=None):
        if gamma is None: gamma = self.gamma
        b = self.bessel(self.msx, gamma*self.rc)
        bsel = abs(b)>1e-12
        
        alphac = zeros_like(self.mdata)
        alphac[bsel] = self.mdata[bsel]/b[bsel]
        return alphac
    alphac = property(calc_alphac)

    def __call__(self, rext, fourier=True):
        "Return external fields calculated at r radial coordinates"
        ms = self.msx[...,newaxis]
        ext = self.alphac[...,newaxis]*self.bessel(ms,self.gamma*rext)
        ext = rollaxis(ext,2)
        
        if fourier:
            return ext
        else:
            extp = fft.ifft(ext, axis=-1)
            return extp*exp(1j*self.m0*self.phi)

    def vector(self, rext, fourier=True):
        "Return r,phi,z vectors reconstructed from +/- coefficients, factor=beta"
        msp,msm = self.msx
        ap, am = self.alphac

        hr = 0.5*(ap*self.bessel(msp, self.gamma*rext) + am*self.bessel(msm, self.gamma*rext))
        hp = -0.5j*(ap*self.bessel(msp, self.gamma*rext) - am*self.bessel(msm, self.gamma*rext))
        hz = 1j*self.gamma/self.factor/2*(ap-am)*self.bessel(self.ms, self.gamma*rext)
        
        return asarray([hr,hp,hz])

class AnalyticExtensionLeft(AnalyticExtension):
    def set_extension_mode(self, mode):
        if self.interior:
            gamma = mode.gamma_int
            self.factor = mode.k0*mode.interior_index**2
        else:
            gamma = mode.gamma_ext
            self.factor = mode.k0*mode.exterior_index**2

        mdata =mode.get_left(field=True, swap=mode.reverse_left)
        return self.set_extension(mdata, mode.m0, mode.pv, mode.coord, gamma)

    def vector(self, rext, fourier=True):
        "Return r,phi,z vectors reconstructed from +/- coefficients, factor=k n^2"
        msp, msm = self.msx
        ap, am = self.alphac

        er = 0.5*(ap*self.bessel(msp, self.gamma*rext) + am*self.bessel(msm, self.gamma*rext))
        ep = -0.5j*(ap*self.bessel(msp, self.gamma*rext) - am*self.bessel(msm, self.gamma*rext))
        ez = 1j*self.gamma/self.factor/2*(ap-am)*self.bessel(self.ms, self.gamma*rext)
        
        return asarray([er,ep,ez])


