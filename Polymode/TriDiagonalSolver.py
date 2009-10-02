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
TriDiagonalSolver.py
===========

Implements a dense block-tridiagonal solver

"""

from __future__ import division
import logging, datetime
from numpy import inf, dtype, array, zeros, append, conj

from . import Material, Waveguide, Equation, Modes
from .Solver import Solve
from .mathlink import blockarray,eigensolver,timer
from .difflounge import finitedifference, boundary

#Import eigensolver
use_arpack=True
if use_arpack:
    #ARPACK interface - this keeps changing in scipy
    from scipy.sparse.linalg.eigen.arpack import eigen
    from scipy.sparse.linalg.interface import LinearOperator
else:
    from .mathlink.eigensolver import eigs

#Try importing external block LU solver
try:
    from .mathlink import ublocklu
except ImportError:
    logging.debug("Warning: Couldn't find C++ block lu solver, using python version")
    ublocklu = None

#**************************************************************************************
# Linear operator routines for eigensolver
#**************************************************************************************
LinearOperator = object

class ShiftInvertPyBlock(LinearOperator):
    def __init__(self, overwrite=False):
        self.lu = blockarray.TriBlockLU(overwrite=overwrite)
        self.overwrite = overwrite
        self.transpose = False

    def set_shift(self, Ain, shift=0.0):
        self.lu(Ain, shift)
        self.shape = (Ain.shape[0],)*2
        self.shift = shift
        self.dtype = Ain.dtype

    def eigenvalue_transform(self, hvals):
        return 1./hvals + self.shift

    def update(self, Aupdate, uprows=1):
        self.lu.update(Aupdate, uprows=uprows)
        
    def matvec(self, x):
        y = self.lu.solve(x)
        return y

    def rmatvec(self, x):
        y = conj(self.lu.solve_transpose(conj(x)))
        return y

#**************************************************************************************

class ShiftInvertCBlock(LinearOperator):
    def __init__(self, overwrite=False):
        self.lu = ublocklu.cblocklu()
        self.lumatrix = None
        self.overwrite = overwrite

    def set_shift(self, Ain, shift=0.0):
        "Setup LU factorization A-sI = LU with a shift s"
        self.yshape = (Ain.mshape[0], Ain.blockshape[0])
        self.shift = shift
        self.shape = (Ain.shape[0], Ain.shape[0])
        self.dtype = Ain.dtype

        if self.overwrite:
            self.lu(Ain, Ain.blockshape[0], shift)
        else:
            self.lumatrix = None
            self.lumatrix = Ain.copy()
            self.lu(self.lumatrix, Ain.blockshape[0], shift)

    def eigenvalue_transform(self, hvals):
        return 1./hvals + self.shift

    def update(self, Aupdate, uprows=1):
        "Update LU factorization A-sI = LU using last uprows of Aupdate"
        self.lu.update(Aupdate, uprows)

    def matvec(self, x):
        y=x.copy()
        self.lu.solve(y)
        return y

    def rmatvec(self, x):
        y=conj(x)
        self.lu.solve_transpose(y)
        return conj(y)

if ublocklu is None:
    ShiftInvertBlock = ShiftInvertPyBlock
else:
    ShiftInvertBlock = ShiftInvertCBlock

#**************************************************************************************

class TriDiBlockSolve(Solve):
    '''
    Main solve class for VWE/SWE. Construct a solver object with
        Solve(shape_in, wg, store=True, compress_to_size=None, discard_vectors=False,
                 mode_calculations=False, label={}, dtype=complex128):
        Where:
            shape (Nr,Naz): The number of radial and Azimuthal nodes to use
            wg: The waveguide
            compress_to_size: size to compress to or None to not store modes
    '''
    overwrite = False       #Overwrite matrix with LU decomp
    stop_test = None        #Give a function that enables a premature stop

    def setup(self, tolerance=1e-10, bandwidth=3, bcwidth=3, xbc=0, \
                fast_convolve=0, eqtype='vector', add_unconverged=False, force_electric=True):
        '''
        Create new Vector Wave Equation with custom parameters:
        
        bandwidth: Finite difference matrix bandwidth
        bcwidth: Width of stencil used for the boundary conditions
                 (must - currently - be less than matrix bandwidth)
        xbc: location of expansion for boundary condition stencil 
        eqtype: either vector or scalar
        '''
        coord = self.wg.get_coord(self.base_shape, border=0)

        #Boundary conditions are Neumann for m+m0==0 and Dirichlet otherwise
        if coord.rmin==0:
            bcl = boundary.BC_Switch(coord.rv, bandwidth=bcwidth, xbc=xbc, dtype=self.dtype)
        else:
            bcl = boundary.BC_Mixed(coord.rv, bandwidth=bcwidth, xbc=xbc, dtype=self.dtype)
        bcr = boundary.BC_Mixed(coord.rv, bandwidth=bcwidth, xbc=xbc, dtype=self.dtype)
        
        #Setup default finite differences
        fd_diff = finitedifference.DifferenceMatrix(2, bandwidth=bandwidth, \
                    X=coord.rv, bcr=bcr, bcl=bcl, dtype=self.dtype)
        fd_diff_ip = finitedifference.DifferenceMatrix(2, bandwidth=bandwidth, \
                    X=coord.rv, bcr=bcr, bcl=bcl, dtype=self.dtype)
        fd_jac = finitedifference.JacobianMatrix(2, bandwidth=bandwidth, \
                    X=coord.rv, bcr=bcr, bcl=bcl, dtype=self.dtype)
        
        if eqtype=='vector':
            eqkwargs = {'fast_convolve':fast_convolve, 'dtype':dtype}
            self.equation = Equation.VectorWaveEquation(fd_diff, **eqkwargs)
            self.equationip = Equation.VectorWaveEquation(fd_diff_ip, **eqkwargs)
            self.jacobian = Equation.VectorWaveJacobian(fd_jac, **eqkwargs)
            self.mode_class = Modes.VectorMode
        
        elif eqtype=='scalar':
            raise NotImplementedError, "Scalar Wave Equation not implemented"
            #self.equation = Equation.Scalar_Wave_Equation(fd_diff)
            #self.jacobian = Equation.Scalar_Wave_Jacobian(fd_jac)
            #self.mode_class = Modes.ScalarMode
        else:
            logging.error("Equation type not recognised")

        #Misc solver parameters
        self.tolerance = tolerance

        #Save unconverged modes if true, else discard them
        self.add_if_unconverged = add_unconverged
        self.force_electric_calculation = force_electric

    # +-----------------------------------------------------------------------+
    # | Pickling marshalling functions
    # | We don't pickle the matrix data to give compact storage
    # +-----------------------------------------------------------------------+

    def __getstate__(self):
        "Pickle all needed data, ignore cached data"
        state = self.__dict__.copy()
        ignore_list = ["matrix", "si"]
        for ignore in ignore_list:
            if ignore in state:
                state[ignore] = None
        return state
    
    def __setstate__(self,state):
        "Restore pickled data"
        self.__dict__.update(state)

    # +-----------------------------------------------------------------------+
    # | Matrix creation functions
    # | These are unnessesary for matrix free methods ..
    # +-----------------------------------------------------------------------+

    def generate(self, leftbc=0, rightbc=0):
        '''
        Generate VWE Matrix
        leftbc, rightbc: only update the specified boundary conditions
        '''
        Nr, Naz = self.base_shape
        bw = self.equation.diff.bandwidth
        pmax = self.equation.pmax
        M = self.matrix.blockview()
        
        #Extents over which the boundary condition changes the lines
        bcstart,bcend = self.equation.diff.bc_extents()
        
        #Only update or full construction
        if not leftbc and not rightbc:
            rows = range(Nr)
        else:
            rows=array([])
            if leftbc:
                rows = append(rows, range(bcstart))
            if rightbc:
                rows = append(rows, range(Nr-bcend, Nr))
        
        #Block constrcution of matrix
        vec_row = zeros((1,Naz,pmax), dtype=self.dtype)
        for ii in rows:
            blockstart = max(bw//2 - ii,0)
            for kk in range(pmax*Naz):
                vec_row.flat[kk] = 1
                y = self.equation.construct(vec_row, ii)
                M[ii,blockstart:blockstart+y.shape[0],:,kk] = \
                    y.reshape((y.shape[0], 2*Naz))
                vec_row.flat[kk] = 0

    def create(self, generate=True):
        '''
        Create the initial equation matrix
        '''
        t = timer.timer(); t.start()
    
        Nr, Naz = self.base_shape
        bw = self.equation.diff.bandwidth
        blockshape = (self.equation.pmax*Naz,)*2
        self.matrix = blockarray.BlockArray((Nr,bw), blockshape=blockshape, dtype=self.dtype)

        #Final setup of the equation objects
        self.equation.setup(self.base_shape,self.wg,self.m0,self.wl)
        self.equationip.setup(self.base_shape,self.wg,self.m0,self.wl)
        self.jacobian.setup(self.base_shape,self.wg,self.m0,self.wl)
        
        #This is a nasty hack!
        self.equationip.diff.generate(leftbc=1, rightbc=1, nodefault=1)
        
        if generate: self.generate()
        logging.debug("Matrix generation done in %.4gs" % t.lap())

    def update_lambda(self, ev, update=1):
        "Set the boundary eigenvalue and update those rows affected"
        self.equation.set_lambda(ev)
        self.jacobian.set_lambda(ev)
        self.generate(leftbc=update, rightbc=update)

    def residue(self, x, l=None):
        if l is None: l = x.evalue

        res = resl = 0
        self.equation.set_lambda(l)
        if x.right is not None:
            res = absolute(self.equation.matvec(x.right)-l*x.right).max() \
                /linalg.norm(x.right)

        if x.left is not None:
            resl = absolute(self.equation.rmatvec(x.left)-conj(l)*x.left).max() \
                /linalg.norm(x.left)

        return max(res, resl)

    def get_data(self):
        return self.modes

    def clear_data(self):
        self.modes = []

    def clean_up(self):
        self.modes = []
        self.matrix = None
        self.si = None

class FixedPointSolver(TriDiBlockSolve):
    """
    Find modes with the real part of the effective index in the given range
    Run with the following parameters:
    m0, wavelength, neffbracket, number=10
    """
    def __init__(self, *args, **kwargs):
        Solve.__init__(self, *args, **kwargs)

    def setup(self, **kwargs):
        tolerance = 1e-8
        iterations = 10
        
        if 'tolerance' in kwargs:
            tolerance = kwargs['tolerance']
        if 'iterations' in kwargs:
            iterations = kwargs.pop('iterations')

        #Iterative parameters
        self.abc_iterations = iterations
        self.abc_convergence = tolerance
        Solve.setup(self, **kwargs)

    def calculate_one(self, evsearch, searchnum):
        """
        Calculate single mode using fixed point iterations
        calculate_one(evsearch, searchnum)
        """
        m0 = self.m0; k0 = self.k0
        wl = self.wl
        
        bcstart, bcend = self.equation.diff.bc_extents()
        coord = self.wg.get_coord(self.base_shape, border=1)

        #Create new mode
        mode = self.mode_class(coord=coord, symmetry=self.wg.symmetry, m0=m0, wl=wl, evalue=evsearch, wg=self.wg)
        mode.right = ones(mode.shape, dtype=self.dtype)
        mode.discard_vectors = self.discard_vectors

        residue = inf; numit=0; evtol=1e-10
        while residue>self.abc_convergence and numit<self.abc_iterations:
            numit+=1
            
            #Update equation with new BCs and update LU decomp
            self.update_lambda(mode.evalue)
            if self.equation.coord.rmin==0:
                self.si.update(self.matrix, uprows=bcend)
            else:
                self.si.set_shift(self.matrix, complex(mode.evalue))

            #Solve linear eigenproblem
            if use_arpack:
                evals, revecs = eigen(self.si, k=searchnum, \
                    which='LM', return_eigenvectors=True)
                evals = self.si.eigenvalue_transform(evals)
            else:
                evals, revecs = eigs(self.si, searchnum, tol=evtol)
                revecs = revecs.T
            
            #Locate closest eigenvalue
            itrack = absolute(evals-mode.evalue).argmin()
            mode.evalue = evals[itrack]
            mode.right = asarray(revecs)[:,itrack]
            
            #Residue and convergence
            evconverge = abs(mode.evalue - evals[itrack])
            residue = mode.residue = self.residue(mode)
            mode.convergence += [ residue ]
            mode.track += [ mode.evalue ]
            
            logging.debug( "[%d] neff: %s conv: %.3g, res:%.3g" % (numit, mode.neff, evconverge, mode.residue) )

        mode.iterations = numit

        #Discard vectors at this stage, to free up memory if we can
        mode.compress()
        
        return mode, evals


class CenterSolver(FixedPointSolver):
    """
    Find modes with the real part of the effective index in the given range
    Note overwrite must be false
    Run with the following parameters:
    m0, wl, neffbracket, number=10
    """
    #Unlike all other solvers overwrite should _default_ to false
    #as the creation of the matrix currently takes so long
    overwrite = False
    def initialize(self, m0, wl, neffrange=None, number=10):
        """
        Setup a solver run with paramters:
        m0, wl, neffs, searchnumber = 1
        """
        FixedPointSolver.initialize(self, m0, wl, neffrange, number)
        self.ignore_outside_interval = size(self.bracket)>1
        
        #Solver specific parameters
        self.searchnumber = 3
        self.numbersolved = 0
    
    def isfinished(self):
        "Check if the solver has found enough modes or is outside the bracket"
        if self.is_finalized: return True

        number_test = (len(self.modes)>0)
        #custom_test = self.user_stop_test(self.modes)
        return number_test
    
    def calculate(self, number=inf):
        m0=self.m0; k0=self.k0
        number = min(number, self.totalnumber)
        
        logging.info( "Center solver. Finding modes, m0=%d, wl=%.3g" % (m0, self.wl) )
        
        #Create new matrix if there is no existing one
        if self.si is None:
            self.si = ShiftInvertBlock(overwrite=self.overwrite)
            self.create()
        
        #Time mode solve
        tick = timer.timer()
        tick.start()

        #Array to store solved modes
        evrange = array(self.bracket)**2*self.k0**2
        
        #Update equation with new BCs and update LU decomp
        self.update_lambda(self.evapprox, not self.overwrite)
        self.si.set_shift(self.matrix, self.evapprox+0j)

        #Find modes of linearized system
        if use_arpack:
            cev = eigen(self.si, k=self.totalnumber, which='LM', return_eigenvectors=False)
            cev = self.si.eigenvalue_transform(cev)
        else:
            cev, cvecs = eigs(self.si, self.totalnumber, tol=1e-6)
        
        #Filter cev within range
        if self.ignore_outside_interval:
            cev = filter(lambda x: min(evrange)<=real(x)<=max(evrange), cev)
        
        #Refine modes
        for ii in range(len(cev)):
            evsearch = cev[ii]
            
            #If the matrix is overwritten we must recreate it
            if self.overwrite: self.generate()
            self.si.set_shift(self.matrix, evsearch+0j)

            #Calculate next mode
            mode, evals = self.calculate_one(evsearch, self.searchnumber)

            #Add mode to list
            if (mode.residue<self.abc_convergence) or self.add_if_unconverged:
                self.modes += [ mode ]

            avtime = tick.lap()/(ii+1)
            logging.info( "Mode #%d [%d/%.3gs], neff=%s, res: %.2e" % \
                (self.numbersolved, mode.iterations, avtime, mode.neff, mode.residue) )

        self.add_time(tick.lap())
            
        #Clean up if calculation is finished!
        if self.isfinished(): self.finalize()

        return self.modes

DefaultSolver = CenterSolver

