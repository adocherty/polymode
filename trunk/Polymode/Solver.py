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
Solver.py
===========

Main solve class for Polymode

Utility functions:

batch_file_save(solvers, filename=None)
 - Save solvers in batch file for later solution

batch_file_load(filename=None)
 - Load solvers in from a batch file

batch_file_load_modes(filename=None)
 - Load modes directly from batch file

batch_solve(solvers, filename=None)
 - Solve problems in list solvers saving periodically to the specified file if given.

batch_continue(filename)
 - Continue an aborted batch_solve with a batch file

export_to_nader(solver, prefix="")
 - Export waveguide.in and input_paramters.in to be read by Nader's solver
"""
from __future__ import division
import logging

import numpy as np
from numpy import *

from . import Material, Waveguide, Equation, Modes

# Cached random functions
class CachedRandom(object):
    """
    Create
    """
    def __init__(self):
        self.cache = []
        self.index = 0

    def reset(self):
        self.index = 0

    def __call__(self, shape):
        from scipy import random
        if self.index>=len(self.cache):
            self.cache += [random.random(shape)]
        x = self.cache[self.index]
        self.index += 1
        return x

#**************************************************************************************

class Solve(object):
    '''
    Main solve class for VWE/SWE. Construct a solver object with
        wg:             The waveguide
        Nres:       The resolution of the calculation grid
        store:      Store mode field information if true
        label:          Dict of information to associate with the solver and mode

        compress_to_size:  Size to compress to or None to not store modes
        mode_calculations: Calculate mode information before discarding fields
    '''
    def __init__(self, wg, Nres=None, store=True, compress_to_size=None,
                 mode_calculations=False, label={}):
        self.wg = wg
        self.base_shape = Nres
        
        #General solver options - All solvers should support these paramters
        self.store_mode_properties = mode_calculations
        self.store_fields = store
        self.compress_to_size = compress_to_size
        self.force_electric_calculation = False
        
        self.dependancies = []              #Don't run this solver until these are true
        self.label = label                          #Custom label to identify the solver
        self.dtype = complex128
        self.modes = []

        #Setup equation with default parameters, can call it to customize
        #Solver specific paramters
        self.setup()

    def setup(self):
        #Solver specific paramteres
        pass

    def add_dependancy(self, depends):
        "Add solver to dependancy list"
        depends = atleast_1d(depends)
        
        for d in depends:
            if hasattr(d,'id'):
                self.dependancies.append(d.id)
            elif 0<int(d)<len(Solve.ids):
                self.dependancies.append(Solve.ids[d])
            else:
                raise LookupError , "Dependancy not recognised, should be a solver"

    # +-----------------------------------------------------------------------+
    # | General solver functions .. may be overloaded
    # +-----------------------------------------------------------------------+

    def plot(self):
        "Plot the effective indices of the found modes"
        import pylab as p_
        col_red = array([0.8,0,0.2])
        col_blue = array([0.2,0,0.8])
        
        neffs = [md.neff for md in self.modes]
        spurious = array([md.guess_spurious() for md in self.modes])
        nconverged = array([md.residue>self.tolerance for md in self.modes])
        colors = col_red*spurious[:,newaxis] + col_blue*nconverged[:,newaxis]
        
        p_.scatter(real(neffs), imag(neffs), s=5, c=colors, marker='o')

    ##
    ## General modal functions
    ##
    def calculate_mode_properties(self, mode):
        #Use rough sampling
        cprop = coord.new(Nshape=(50, 5*mode.symmetry), symmetry=1)
        #cprop = coord.new(Nshape=(100, 11), symmetry=mode.symmetry)

        calc = {}
        calc['effective_area'] = mode.effective_area(coord=cprop)
        calc['spot_size'] = mode.spot_size(coord=cprop)
        calc['group_velocity'] = mode.group_velocity(self.wg, coord=cprop)

        mode.calculated_info = calc

    def guess_spurious_mode(self, mode, cutoff=5.0):
        '''
        Guess if this is a real or spurious mode based on the mode
        intensity distribution
        '''
        #The outermost object in the waveguide, guess if not given
        router = 0.95*self.wg.get_rmax(0)
        c = mode.coord

        #If the coord object doesn't have a rv member or
        #mode doesn't have field information this will fail
        try:
            #RMS of magnetic intensity over azimuthal direction
            hr,ha,hz = mode.magnetic_field()
            pprofile = mean(abs(hr)**2+abs(ha)**2+abs(ha)**2,axis=1)
            pfrac = mean(pprofile[c.rv>router])/mean(pprofile[c.rv<router])
        except:
            pfrac=0
            
        mode.is_spurious = pfrac>cutoff
        return mode.is_spurious

    def residue(self, x, l=None):
        pass

    ##
    ## Interface to MPI solver commands
    ##
    def get_data(self):
        return self.modes

    def clear_data(self):
        self.modes = []

    def clean_up(self):
        self.modes = []

    def calculate(self, number=inf):
        pass

    def __call__(self, *args, **kwargs):
        '''
        Solve the constructed problem with
        m0: the waveguide symmetry index
        wl: wavelength
        neffrange: the upper and lower real effective indices of the search range
        nefflist: Find modes near these effective indices
        modelist: Find modes near these modes
        totalnumber: total number of modes to find
        '''
        self.initialize(*args, **kwargs)
        self.calculate()
        self.finalize()
        return self.modes
        
    def initialize(self, wl, m0=0, neffrange=None, nefflist=None, modelist=None, number=1):
        '''
        Setup the solver with pre-calculation parameters with:
        wl: wavelength
        m0: the waveguide symmetry index
        neffrange: the upper and lower real effective indices of the search range
        nefflist: Find modes near these effective indices
        modelist: Find modes near these modes
        totalnumber: total number of modes to find
        '''
        self.m0 = m0
        self.wl = wl
        self.k0 = 2*pi/wl

        #Set number and neffrange depending on the case
        self.numbercalculated = 0
        if nefflist is not None:
            self.bracket = 0,inf
            self.totalnumber = len(nefflist)
        elif modelist is not None:
            self.bracket = 0,inf
            self.totalnumber = len(modelist)
        else:
            #Calculate range from core index if not given
            self.bracket = self.wg.index_range(wl)
            self.totalnumber = number

        #Or manual setting
        if neffrange is not None:
            if iterable(neffrange):
                self.bracket = neffrange
            else:
                self.bracket = (self.wg.index_range(wl)[0], neffrange)

        self.clean_up()

        #Mode/neff lists if any
        self.nefflist = nefflist
        self.modelist = modelist
        self.is_finalized = False

    def estimate_complete_fraction(self):
        "Return a number between 0 (started) and 1 (finished)"
        return float(len(self.modes))/self.totalnumber

    def finalize(self):
        """
        Finalize the modes after the solver has finished.
        
        Including
         - Clean up temporary objects
         - Delete or compress mode vectors is required
         - Remove debug information if not in debugging mode
        """
        logging.info("Finalizing calculated modes")
        for ii,mode in enumerate(self.modes):
            #Label the mode
            mode.label = self.label

            #Store some quantities we may want if we throw away the vector info
            if self.store_mode_properties:
                self.calculate_mode_properties(mode)

            #Update spurious indicator
            self.guess_spurious_mode(mode)

            #Remove calculated EF if forced
            if self.store_fields:
                mode.store_calculated_electric_field(wg=self.wg, force=self.force_electric_calculation)

                if self.compress_to_size is not None:
                    mode.compress(self.wg, self.compress_to_size)

                #Add extension for behaviour outside the computational domain
                mode.normalize(wg=self.wg)
            else:
                mode.discard_fields()

        #Sort modes
        self.modes.sort(reverse=True)
        self.is_finalized = True


class WavelengthTrack(Solve):
    '''
    Track modes over a wavelength range with adaptive step size
    '''
    def __init__(self, solver, track_range=None):
        self.solver = solver
        self.track_range = None
        self.ga_target = 1e-3

        Solve.__init__(self, solver.wg, compress_to_size=solver.compress_to_size)
        
    def initialize(self, wl_range, m0=0, **kwargs):
        self.wl_range = wl_range
        
        self.m0 = m0
        self.solver.initialize(wl_range[0], m0, **kwargs)

    def calculate(self, number=inf):
        import pylab as pl
        solver = self.solver

        #Setup wavelength range
        wl_start, wl_stop = self.wl_range

        #Starting step size
        dwl = (wl_stop-wl_start)/100.0
        
        #Tolerances for adaptive step sizes
        dwl_minimum = dwl/10
        dwl_maximum = 5*dwl

        ga_target = self.ga_target
        ga_minimum = ga_target/10
        ga_maximum = ga_target*10

        #Find start modes to track
        modes = self.solver.calculate()
        
        #Tracking modes
        Nm = len(modes)
        dneffdwl = zeros(Nm, complex_)
        modes_track = [m.copy() for m in modes]

        num_eval_backtrack = num_eval = 0
        do_update=True

        wl = wl_start
        self.modes = list(modes)
        while wl<wl_stop:
            #Update wavelength
            wl += dwl
            logging.info("")
            logging.info("WL %.6g, step size: %.4g" % (wl,dwl))

            #Find new modes
            modes_last = modes
            self.solver.initialize(wl, self.m0, modelist=modes_track)
            modes = self.solver.calculate()
            num_eval +=1

            if 0:
                print "Neff guess accuracy:", abs(modes[0].neff-neff_guess)/abs(neff_guess)

                m1 = modes[0]
                solver.equation.set_lambda(m1.evalue)
                M1x = solver.equation.matvec(m1.right) - m1.evalue*m1.right
            
                solver.jacobian.setup(solver.base_shape,solver.wg,self.m0,wl+dwl)
                solver.jacobian.set_lambda(m1.evalue)
                M0px = solver.jacobian.matvec(m1.right) - m1.right

                dmu = -dot(conj(m1.left), M1x)/dot(conj(m1.left), M0px)
                neff_guess = sqrt(m1.evalue+dmu)/(2*pi/m1.wl)
            
            Nm_current = len(modes)
            if Nm_current==0:  #Jump to next point and try and find modes there 
                continue
            
            elif Nm_current<Nm:  #Find a replacement mode?
                print "Lost %d modes" % (Nm- Nm_current)

            elif Nm_current>len(modes_last):
                logging.warning("Found more modes than requested!")
                modes = modes[:len(modes_last)]
               
            #Calculate mode differences
            remove_modes = []
            dneffdwl_last = dneffdwl
            dneffdwl = zeros(Nm_current, complex_)
            ga_max = 0; ga_min = inf
            for ii in range(Nm_current):
                neff = modes[ii].neff

                #Calculate dispersion
                dneffdwl[ii] = (neff - modes_last[ii].neff)/dwl

                #Guess accuracy
                ga = abs(modes_track[ii].neff-neff)/abs(modes_track[ii].neff)
                ga_max=max(ga_max,ga); ga_min=min(ga_min,ga)

                #Have the modes left the tracked range?
                if self.track_range is not None and (neff<min(track_range) or neff>max(track_range)):
                    logging.warning("Mode has left tracked neff range")
                    remove_modes.append(ii)

            #Adaptive guess for next dwl
            accept = True
            if wl>wl_start+dwl:
                if ga_max>0:
                    dwl_target = dwl*(ga_target/ga_max)**(0.5)

                if (ga_max>ga_maximum) and (dwl>dwl_minimum):
                    logging.info("Eigenvalue change to large. Backtracking")
                    accept = False
                    dwl_target = min(dwl_target,dwl*0.5)

                dwl = dwl_target

            #Guess next neff
            if accept:
                dneffdwl_last = dneffdwl
                self.modes += modes
            #Backtrack!
            else:
                wl -= dwl
                modes = modes_last
                dneffdwl = dneffdwl_last
                num_eval_backtrack +=1
            
            #Truncate modes_track otherwise modes can be larger than modes_last
            modes_track = [m.copy() for m in modes]
            
            #Update neff for modes_track
            for ii in range(len(modes)):
                modes_track[ii].neff = (modes[ii].neff + dneffdwl[ii]*dwl)

            logging.info("Dispersion: %s " % dneffdwl)
            logging.info("Guess accuracy: %0.4g -> %0.4g" % (ga_max, ga_min))

        logging.info("Total points: %d, number of backtracks: %d" % (num_eval, num_eval_backtrack))
        return self.modes
        
    def update_eigenvector(self,m1,m2):
            #Calculate perturbation
#           eps = 1e-3
#           solver.equation.setup(solver.base_shape,solver.wg,m0,m1.wavelength+eps)
#           solver.equation.set_lambda(m1.evalue)
#           M1xp = solver.equation.matvec(m1.right)

#           solver.equation.setup(solver.base_shape,solver.wg,m0,m1.wavelength-eps)
#           solver.equation.set_lambda(m1.evalue)
#           M1xm = solver.equation.matvec(m1.right)
#           M1x = (M1xp-M1xm)/(2*eps)

            solver.equation.setup(solver.base_shape,solver.wg,m0,m2.wavelength)
            solver.equation.set_lambda(m1.evalue)
            M1x = solver.equation.matvec(m1.right) - m1.evalue*m1.right
            
            solver.jacobian.setup(solver.base_shape,solver.wg,m0,m1.wavelength)
            solver.jacobian.set_lambda(m1.evalue)
            M0px = solver.jacobian.matvec(m1.right) - m1.right

            dmu = -dot(conj(m1.left), M1x)/dot(conj(m1.left), M0px)
            dneffc1 = (m1.neff**2/m1.wavelength+0.5*dmu/m1.k0)/m1.neff
            dneffc = sqrt(m1.evalue+dmu)/m2.k0 - m1.neff
            print "dneff(1)", dneffc1
            print "dneff(2)", dneffc
            print
            
            neff_guess += [sqrt(m1.evalue+dmu)/m2.k0]
            
            #Find correction to eigenvector
            mu2 = m2.evalue+0*dmu
            Mx1 = -(M0px*dmu/delta + M1x)
            
            #Approx:
            if not hasattr(solver, 'matrix'):
                Nr, Naz = solver.base_shape
                bw = solver.equation.diff.bandwidth
                blockshape = (solver.equation.pmax*Naz,)*2
                solver.matrix = blockarray.BlockArray((Nr,bw), blockshape=blockshape, dtype=complex_)
                si = Solver.ShiftInvertBlock(overwrite=False)
            
            solver.generate()
            si.set_shift(solver.matrix, complex(m1.evalue))
            x1 = si.matvec(Mx1)
            
            y = m1.right + delta*x1

            solver.equation.set_lambda(m2.evalue)
            print "Diff1", linalg.norm(solver.equation(y)-m2.evalue*y)
            print "Diff2", linalg.norm(solver.equation(m1.right)-m2.evalue*m1.right)
            
            #plot([(m1.wavelength+m2.wavelength)/2], [dneffm], 'rs')
            #plot([m1.wavelength], [dneffc], 'bo')


class WavelengthScan(Solve):
    '''
    Scan over a wavelength range and plot a condition number for the modal
    eigenvalue problem. The exact nature of this condition number depends
    upon the nature of the algorithm in the supplied solver
    '''
    def __init__(self, solver, Nscan=(20,100)):
        self.solver = solver
        self.Nscan = Nscan

        #The condition number scan is stored here
        self.Cscan = np.zeros(self.Nscan, dtype=float)
        self.neffscan = np.zeros(self.Nscan, dtype=float)
        
        Solve.__init__(self, solver.wg, compress_to_size=solver.compress_to_size)
    
    def initialize(self, wl_range, m0=0, **kwargs):
        self.wl_range = wl_range
        
        self.m0 = m0
        self.solver.initialize(wl_range[0], m0, **kwargs)

        if 'neffrange' in kwargs:
            self.neffrange = kwargs['neffrange']
        else:
            self.neffrange = None

    def calculate(self, number=inf):
        import pylab as pl
        solver = self.solver

        #Setup wavelengths
        dwl = (self.wl_range[1]-self.wl_range[0])/self.Nscan[0]
        wls = np.arange(self.wl_range[0], self.wl_range[1], dwl)

        for ii,wl in enumerate(wls):
           #Update wavelength
            self.solver.initialize(wl, self.m0)

            #Range to scan
            if self.neffrange is None:
                neffrange=self.wg.index_range(wl)
            else:
                neffrange=self.neffrange

            dneff = (neffrange[1]-neffrange[0])/self.Nscan[1]
            neffs = np.arange(neffrange[0], neffrange[1], dneff)
            
            #Scan over beta range
            self.Cscan[ii] = self.solver.condition(neffs*self.solver.k0)
            self.neffscan[ii] = neffs

        return self.Cscan
        
    def plot(self, style={}):
        import pylab as pl

        dwl = (self.wl_range[1]-self.wl_range[0])/self.Nscan[0]
        wls = np.arange(self.wl_range[0], self.wl_range[1], dwl)
        wlscan = wls[:,newaxis] + 0*self.neffscan
        
        pl.pcolor(wlscan, self.neffscan, np.log10(self.Cscan), **style)

        if 0:
            pl.plot(betascan/self.solver.k0, self.Cscan[ii])



def batch_file_save(solvers, filename=None):
    "Save solvers in batch file for later solution"
    from cPickle import dump
    try:
        dump(solvers, open(filename,'wb'))
    except IOError:
        logging.error("Failed to save solvers to file %s" % filename)

def batch_file_load(filename=None):
    "Load solvers in from a batch file"
    from cPickle import load
    try:
        solvers = load(open(filename,'rb'))
    except IOError:
        solvers = []
        logging.error("Failed to load batch solver file %s" % filename)
    return solvers

def batch_file_load_modes(filename=None, return_wg=False):
    "Load modes from batch file"
    solvers = batch_file_load(filename)

    #Add modes to list
    modes = []
    wgs = []
    for solver in solvers:
        modes += solver.get_data()          #This must return a list!
        wgs.append( solver.wg )

    #Return waveguides if requested
    if return_wg:
        return modes, wgs
    else:
        return modes

def batch_solve(solvers, filename=None):
    """
    Solve problems in list solvers saving periodically
    to the specified file if given.
    The batch solve can be continued if interrupted
    with the function `batch_continue(filename)`.
    """
    from cPickle import dump
    for solver in solvers:
        #Resume calculation if not finished
        if not solver.isfinished():
            solver.calculate()

        #Save solver queue
        if filename is not None:
            dump(solvers, open(filename,'wb'))

    modes = []
    for solver in solvers:
        modes += solver.get_data()

    return modes

def batch_continue(filename):
    """
    Continue an aborted batch_solve with a batch file
    """
    solvers = batch_file_load(filename)
    return batch_solve(solvers, filename)

