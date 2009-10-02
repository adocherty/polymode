# _*_ coding=utf-8 _*_
#---------------------------------------------------------------------------------
#Copyright © 2009 Andrew Docherty

#This program is part of Polymode.
#Polymode is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.
#---------------------------------------------------------------------------------

"""from scipy import optimize


Nonlinear Eigenvalue Problem solvers

ToDo:
* Change Mdtn implementation - not very robust at the moment.

"""
from __future__ import division
import logging, datetime

from numpy import *
from numpy.lib.scimath import sqrt
from scipy import optimize

from . import Material, Waveguide, Equation, Modes
from .TriDiagonalSolver import Solve,TriDiBlockSolve,ShiftInvertBlock
from .mathlink import blockarray,eigensolver,timer, utf8out
from .mathlink.misc import format_complex


class NlepResidualIteration(TriDiBlockSolve):
    maxiter = 100
    reshift = 5
    
    # NLEP solver tweaking parameters
    rayleigh_iterations = 4
    
    eo_max_iter = 1
    eo_tolerance = 1e-6
    eo_localize = 0
    eo_maximum = 40
    eo_keep = 35

    evapprox_choice = 'none'
    optimize = 0
    use_real_guess = 1
    use_fastip = 0

    def __init__(self, *args, **kwargs):
        Solve.__init__(self, *args, **kwargs)
        self.overwrite=0

    def ip_dot(self, x1, x2, ev1=None, ev2=None):
        return dot(conj(x1.left),x2.right)

    def ip(self, x1, x2, ev1=None, ev2=None):
        if ev1 is None: ev1 = x1.evalue
        if ev2 is None: ev2 = x2.evalue

        if abs(ev1-ev2)<1e-10:
            self.jacobian.set_lambda(ev2)
            v2 = self.jacobian.matvec(x2)-x2.right
            nip = dot(conj(x1.left), v2)
        else:
            self.equation.set_lambda(ev1)
            v = self.equation.matvec(x2.right)
            self.equation.set_lambda(ev2)
            v -= self.equation.matvec(x2.right)
            nip = dot(conj(x1.left),v)/(ev1-ev2) - dot(conj(x1.left),x2.right)
        return nip

#   @timer.time_function(prefix="NL:")
    def ip_fast(self, x1, x2, ev1=None, ev2=None):
        if ev1 is None: ev1 = x1.evalue
        if ev2 is None: ev2 = x2.evalue
        
        a1 = dot(conj(x1.left), x2.right)
        
        #If ev1 and ev2 are almost the same
        if abs(ev1-ev2)<self.tolerance:
            return a1
        
        M1 = self.equation.calc_rightbc_dtn(ev1)
        M2 = self.equation.calc_rightbc_dtn(ev2)
        
        alpha = 2*(1+0.5/self.equation.coord.Nr)/self.equation.coord.dr/(ev1-ev2)
        a2 = sum(conj(x1.shape_vector(x1.left)[-1])*(M1-M2)*x2.shape_vector(x2.right)[-1])
        return a1-alpha*a2

    @timer.time_function(prefix="NL:")
    def rayleigh_quotient(self, x, ev):
        #Rayleigh Quotient is one step of Newton's iteration
        #to solve y* T(l) x = 0 for l
        ii=0; revchange=inf
        while revchange>1e-12 and ii<self.rayleigh_iterations:
            ii+=1
            self.jacobian.set_lambda(ev)
            self.equation.set_lambda(ev)

            Au = self.equation.matvec(x.right)-ev*x.right
            Apu = self.jacobian.matvec(x.right)-x.right
            ev_new = ev - dot(conj(x.left), Au)/dot(conj(x.left), Apu)

            revchange = abs(ev_new-ev)
            ev = ev_new
        return ev

    def linear_rayleigh_quotient(self, x, ev):
        #Rayleigh Quotient is one step of Newton's iteration
        #to solve y*T(l)x/y*x = 0 for l
        self.equation.set_lambda(ev)
        Au = self.equation.matvec(x.right)-ev*x.right
        ev = dot(conj(x.left), Au)/dot(conj(x.left), x.right)
        return ev
    def rayleigh_objective_ri(self, evri, x):
        ev = evri[0]+1j*evri[1]
        self.equation.set_lambda(ev)
        Au = self.equation.matvec(x.right)-ev*x.right
        f = dot(conj(x.left), Au)
        return (f.real, f.imag)
    def rayleigh_jacobian_ri(self, evri, x):
        ev = evri[0]+1j*evri[1]
        self.jacobian.set_lambda(ev)
        Au = self.jacobian.matvec(x.right)-x.right
        f = dot(conj(x.left), Au)
        return array([[f.real, -f.imag],[f.imag, f.real]])
    def rayleigh_objective(self, ev, x):
        self.equation.set_lambda(ev)
        Au = self.equation.matvec(x.right)-ev*x.right
        f = dot(conj(x.left), Au)
        return f
    def minres_objective(self, evri, x):
        ev = evri[0]+1j*evri[1]
        self.equation.set_lambda(ev)
        Tx = self.equation.matvec(x.right)-ev*x.right
        return linalg.norm(Tx)

    def construct_gram_matrix(self, modes, ev=None):
        nmodes = len(modes)
        GL = empty((nmodes,nmodes), dtype=self.dtype)
        GR = empty((nmodes,nmodes), dtype=self.dtype)
        
        self.jacobian.set_lambda(ev)
        self.equation.set_lambda(ev)
        for ii in range(nmodes):
            m1 = modes[ii]
            for jj in range(nmodes):
                m2 = modes[jj]
                Tx = self.equation.matvec(m2.right)-ev*m2.right
                tip = dot(conj(m1.left),Tx)
                GR[ii,jj] = tip/(m1.evalue-ev)
                GL[jj,ii] = tip/(m2.evalue-ev)
        return GL,GR

    def construct_gram_matrix_fast(self, modes, ev=None):
        nmodes = len(modes)
        GL = empty((nmodes,nmodes), dtype=self.dtype)
        GR = empty((nmodes,nmodes), dtype=self.dtype)
        
        M = self.equation.calc_rightbc_dtn(ev)
        alpha = 2*(1+0.5/self.equation.coord.Nr)/self.equation.coord.dr

        for ii in range(nmodes):
            m1 = modes[ii]; Mi = m1.Mdtn
            for jj in range(nmodes):
                m2 = modes[jj]; Mj = m2.Mdtn
                dip = dot(conj(m1.left), m2.right)
                fl = alpha*sum(conj(m1.shape_vector(m1.left)[-1])*(M-Mj)*m2.shape_vector(m2.right)[-1])
                fr= alpha*sum(conj(m1.shape_vector(m1.left)[-1])*(M-Mi)*m2.shape_vector(m2.right)[-1])
                
                GR[ii,jj] = dip+fr/(m1.evalue-ev)
                GL[jj,ii] = dip+fl/(m2.evalue-ev)
        return GL,GR

    def construct_b(self, modes, x, ev):
        nmodes = len(modes)
        bl = zeros((nmodes), dtype=self.dtype)
        br = zeros((nmodes), dtype=self.dtype)
        M = self.equation.calc_rightbc_dtn(ev)
        
        alpha = 2*(1+0.5/self.equation.coord.Nr)/self.equation.coord.dr
        for ii in range(nmodes):
            m1 = modes[ii]
            fl = alpha*sum(conj(x.shape_vector(x.left)[-1])*(M-m1.Mdtn)*m1.shape_vector(m1.right)[-1])
            fr = alpha*sum(conj(m1.shape_vector(m1.left)[-1])*(m1.Mdtn-M)*x.shape_vector(x.right)[-1])
                
            bl[ii] = dot(conj(x.left), m1.right) + fl/(m1.evalue-ev)
            br[ii] = dot(conj(m1.left), x.right) + fr/(ev-m1.evalue)
        return bl,br

    def approximate_orthogonalize(self, omodes, mode, sig=0):
        nmodes = len(omodes)
        if nmodes<1: return (0,0)

        GL,GR = self.construct_gram_matrix_fast(omodes, ev=mode.evalue)
        bl,br = self.construct_b(omodes, mode, mode.evalue)
        
        alphal = conj(linalg.solve(GL,bl))
        alphar = linalg.solve(GR,br)
        return alphal, alphar

    def construct_alpha_mode(self, omodes, mode, alphal, alphar):
        for ii in range(len(omodes)):
            mode.left -= alphal[ii]*omodes[ii].left
            mode.right -= alphar[ii]*omodes[ii].right
        return mode
        
    def orthogonalize_all(self, modes):
        "Orthogonalize modes to be a bi-orthonormal set"
        nmodes = len(modes)
        for ii in range(nmodes):
            m1 = modes[ii]
            for jj in range(0,ii):
                m2 = modes[jj]
                dlr = self.ip_fast(m1, m1, m1.evalue, m2.evalue)
                alpha = self.ip_fast(m1, m2, m1.evalue, m2.evalue)/dlr
                beta = self.ip_fast(m2, m1, m1.evalue, m2.evalue)/dlr
                m1.right -= m2.right*alpha
                m1.left -= m2.left*conj(beta)
            self.normalize(m1)
        return modes

    def normalize(self, mode):
        #Fix norm of right vector to be one
        nalpha=linalg.norm(mode.right)
        ngamma = dot(mode.right,conj(mode.left))
        nbeta = conj(ngamma/nalpha)
        
        #Check for zero vectors
        if ngamma==0: mode.left[0]=mode.right[0]=1; nlr=1

        mode.right *= 1/nalpha
        mode.left *= 1/nbeta

    def get_data(self):
        return self.modes

    def clear_data(self):
        self.modes = []
        self.ssmodes = []

    def iiteration(self, ssmodes, mode):
        nmodes = len(ssmodes)
        evstart = mode.evalue
    
        self.update_lambda(mode.evalue)
        self.si.set_shift(self.matrix, complex(mode.evalue))

        #Linear inverse interations to initialize the eigenvector
        number_initialize = 5
        logging.debug("Initial %d linear iterations" % number_initialize)
        for ii in range(number_initialize):
            mode.right = self.si.matvec(mode.right)
            mode.left = self.si.rmatvec(mode.left)

        stagnation = 0;     locked_on = 0;      sig=0
        niter=0;            res=inf;            over=0
        while niter<self.maxiter and res>self.tolerance:
            #Reshift at start and occasionally after this
            #and once when mode is locked on
            if locked_on==1 or mod(niter,self.reshift)==0:
                if locked_on==1: locked_on = 2
                self.update_lambda(mode.evalue)
                self.si.set_shift(self.matrix, complex(mode.evalue))
            
            #Note: we stop orthogonalization when evchange is less than the tolerance
            #to prevent stagnation due to imprecise orthogonalization (why?)
            eoiter = 0; eoortho = -1
            evchange = inf
            while evchange>self.eo_tolerance and eoiter<self.eo_max_iter:
                eoiter += 1
                evapprox = self.rayleigh_quotient(mode, mode.evalue)
                evchange = abs(evapprox-mode.evalue)

                #Reject large sudden changes in eigenvalue
                if evchange/abs(mode.evalue)>10:
                    logging.debug("Rejecting large eigenvalue change to %s" % evapprox)
                    mode.evalue = mode.evalue - 1e-4
                else:
                    mode.evalue = evapprox

                #Orthogonalize if not locked on
                if not locked_on and len(ssmodes)>0:
                    alphal,alphar = self.approximate_orthogonalize(ssmodes, mode)
                    mode = self.construct_alpha_mode(ssmodes, mode, alphal,alphar)
                    eoortho = max(abs(alphal)) + max(abs(alphar))

                #self.normalize(mode)
            
            #Detect when a mode is converged to the point that orthogonalization
            #should be suspended. This may result in repeated modes, but otherwise
            #it can cause stagnation
            if not locked_on and (niter>self.reshift) and (eoortho>0) and (abs(evchange)<1e-10):
                locked_on = 1

            #Residual inverse iteration updates
            self.equation.set_lambda(mode.evalue)
            du = self.equation.matvec(mode.right)-mode.evalue*mode.right
            dv = self.equation.rmatvec(mode.left)-conj(mode.evalue)*mode.left
            mode.right -= self.si.matvec(du)
            mode.left -= self.si.rmatvec(dv)

            if any(isnan(mode.right)) or any(isnan(mode.left)):
                logging.warning("Nan found in iterate, rejecting")
                return evstart
            
            self.normalize(mode)
            
            #Measure "stagnation"
            Nstag = 4
            if niter>(Nstag+1):
                stagnation = average(array(mode.convergence[-(Nstag-1):]) \
                    /array(mode.convergence[-Nstag:-1]))
            
            #Convergence information
            mode.residue = res = max(abs(du).max(), abs(dv).max())
            mode.track += [mode.evalue]
            mode.convergence += [mode.residue]

            status_str = {0:"", 1:"L", 2:"*"}[locked_on]
            logging.debug( utf8out(u"[%d:%s] %s +- %.2g s:%.2g eo:%.2g ei:%d res:%.2g" \
                % (niter,status_str,mode.evalue,evchange,stagnation,eoortho,eoiter,res)) )
            niter+=1

        return mode.evalue



# ***********************************************************************

class NLResidualSolver(NlepResidualIteration):
    '''
    Extended Residual Inverse Iteration method
    Solves the problem with non-reflective boundary conditions
    for multiple modes using a residual nonlinear eigenvalue solver
    with approximate deflation
    '''
    def isfinished(self):
        "Check if the solver has found enough modes or is outside the bracket"
        if self.is_finalized: return True
        
        number_test = self.numbercalculated>=self.totalnumber
        if len(self.modes)>1:
            #custom_test = self.user_stop_test(self.modes)
            neffs = real([m.neff for m in self.modes])

            if isscalar(self.bracket) or len(self.bracket)<2:
                outside_bracket = 0
            else:
                outside_bracket = sum(neffs>max(self.bracket))
                outside_bracket += sum(neffs<min(self.bracket))

            return number_test or (outside_bracket>(2))
        else:
            return number_test

    def calculate(self, number=inf):
        if self.totalnumber is inf and number is inf:
            find_number_info = 'all'
        else:
            find_number_info = min(self.totalnumber, number)
        
        logging.info( "NLII solver. Finding %s modes, m0=%d, wl=%.3g" \
            % (find_number_info, self.m0, self.wl) )

        #Set this here for now...
        self.ssmodes = []

        #Create matrix and shift invert solver if needed
        if not hasattr(self,'si') or self.si is None:
            self.create()
            self.si = ShiftInvertBlock(overwrite=False)

        #Set the number to the minimum of number and totalnumber
        number = int(min(self.totalnumber, number))
        neffrange = self.bracket
        neffapprox = complex(self.bracket[1]*(1.0-1e-5))

        #Setup initial parameters
        coord = self.equation.coord
        m0 = self.m0
        wl = self.wl

        #Time mode solve
        tick = timer.timer()
        tick.start()

        #Use random cache:
        if hasattr(self,'random'):
            random=self.random
        else:
            from numpy.random import random
        
        #Tag mode to compress vectors to final size
        if type(self.compress_to_size) is tuple:
            coord_compress = self.wg.get_coord(self.compress_to_size, border=1)
        else:
            coord_compress = None

        #Main loop
        ii=0
        isfinished=0
        while not isfinished and ii<number:
            #Choose new neffapprox if in list mode
            if self.nefflist is not None:
                neffapprox = self.nefflist[ii]
            elif self.modelist is not None:
                neffapprox = self.modelist[ii].neff
        
            logging.debug("Finding mode near %s" % neffapprox)
        
            evapprox = (neffapprox*self.k0)**2

            #Add new mode to subspace
            nextmode = self.mode_class(coord=coord, symmetry=self.wg.symmetry, \
                m0=m0, wl=wl, evalue=evapprox, wg=self.wg)
            nextmode.coord_compress = coord_compress
            nextmode.track=[]
            nextmode.residue = inf

            if self.modelist is not None:
                nextmode.right = self.modelist[ii].right
                nextmode.left = self.modelist[ii].left
            else:
                nextmode.right = random(prod(nextmode.shape)).astype(self.dtype)
                nextmode.left = random(prod(nextmode.shape)).astype(self.dtype)
            
            #Iterate on this mode orthogonalizing against previous modes
            self.iiteration(self.ssmodes, nextmode)
            nextmode.solve_time = tick.lapdelta()
            
            #Add dtn coefficiencts for fast inner products
            nextmode.Mdtn = self.equation.calc_rightbc_dtn(nextmode.evalue)

            #Reject mode from ssmodes list if it hasn't converged
            if nextmode.residue>self.tolerance:
                logging.info("Rejecting unconverged mode" )
                nextmode.status = "Rejected"
                if self.add_if_unconverged: self.modes += [nextmode]
            else:
                self.ssmodes += [nextmode]
                self.modes += [nextmode]

            #If we have found more than the maximum modes, resort and shorten the list
            if len(self.ssmodes)>=self.eo_maximum:
                sortinx = argsort(self.ssmodes)
                self.ssmodes = [self.ssmodes[inx] for inx in sortinx[:self.eo_keep]][::-1]
                evapprox = real(self.ssmodes[-1].evalue)
            
            #Next approximate, choose the newmode evalue or the nextmode evalue
            if self.evapprox_choice=='next':
                neffapprox = nextmode.neff
            if self.use_real_guess:
                neffapprox = real(neffapprox)+0j
            
            #Print information about the found mode
            iterations=len(nextmode.track)
            logging.info( "Mode #%s [%d:%.3gs] neff:%s res: %.4g"
                % (ii, iterations, nextmode.solve_time, \
                    format_complex(nextmode.neff), nextmode.residue) )

            ii+=1
            self.numbercalculated+=1

        return self.modes

DefaultSolver = NLResidualSolver

