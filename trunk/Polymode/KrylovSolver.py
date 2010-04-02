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

"""
Nonlinear Eigenvalue Problem solvers for the Polymode library

ToDo:
* Change Mdtn implementation - not very robust at the moment.

"""
from __future__ import division
import logging, datetime

from numpy import *
from numpy.lib.scimath import sqrt
from scipy import linalg as splinalg

from . import Material, Waveguide, Equation, Modes

from .mathlink import blockarray,eigensolver,timer, utf8out
from .Solver import *
from .TriDiagonalSolver import *

from numpy.random import random

class KrylovIteration(TriDiBlockSolve):
    maxiter = 100
    reshift = None

    rayleigh_iterations = 5
    evordering = 1

    def rayleigh_quotient(self, ev, x, y=None):
        #Rayleigh Quotient is one step of Newton's iteration
        #to solve y* T(l) x = 0 for l
        
        #Construct approximate left vector
        if y is None:
            y = zeros_like(x); y[0] = 1
            y = self.si.rmatvec(y)
        
        ii=0; revchange=inf
        while revchange>1e-12 and ii<self.rayleigh_iterations:
            ii+=1
            self.jacobian.set_lambda(ev)
            self.equation.set_lambda(ev)

            Au = self.equation.matvec(x)-ev*x
            Apu = self.jacobian.matvec(x)-x
            ev_new = ev - dot(conj(y), Au)/dot(conj(y), Apu)

            revchange = abs(ev_new-ev)
            ev = ev_new
        return ev

    def approximate_orthogonalize(self, u, v, mu):
        No = self.Northogonalize
        vs = self.rorthogonalize; ws = self.lorthogonalize

        self.equation.set_lambda(mu)
        G = zeros((No,No), complex_)
        bx = zeros((No), complex_)
        by = zeros((No), complex_)

        #Build projected matrix
        Tv = self.equation.matvec(v)
        THu = self.equation.rmatvec(u)
        for ii in range(No):
            Tvii = self.equation.matvec(vs[ii])
            for jj in range(No):
                G[ii,jj] = dot(conj(ws[jj]), Tvii) - mu*dot(conj(ws[jj]), vs[ii])
            
            bx[ii] = dot(conj(ws[ii]), Tv) - mu*dot(conj(ws[ii]), v)
            by[ii] = dot(conj(THu), vs[ii]) - mu*dot(conj(u), vs[ii])
        
        alphar = linalg.solve(G.T,bx)
        alphal = conj(linalg.solve(G,by))
        print alphar, alphal
        for ii in range(self.Northogonalize):
            v -= alphar[ii]*self.rorthogonalize[ii]
            u -= alphal[ii]*self.lorthogonalize[ii]
            
        return u,v

    def construct_alpha_vectors(self, u, v, alphal, alphar):
        for ii in range(self.Northogonalize):
            v -= alphar[ii]*self.rorthogonalize[ii]
            u -= alphal[ii]*self.lorthogonalize[ii]
        return u,v
    
    def get_data(self):
        return self.modes

    def clear_data(self):
        self.modes = []
        self.subspace = []

    def projected_matvec(self, vs, mu, x, prime=False):
        Nv = len(vs)

        vp = zeros(self.vector_shape, complex_)
        Ax = zeros_like(x)
        for ii in range(Nv):
            vp += vs[ii]*x[ii]

        if prime:
            Tv = self.jacobian.matvec(vp) - vp
        else:
            Tv = self.equation.matvec(vp) - mu*vp
        
        for jj in range(Nv):
            Ax[jj] = dot(conj(vs[jj]), Tv)
        return Ax

    def projected_rmatvec(self, vs, mu, y, prime=False):
        Nv = len(vs)

        #Projected vector
        vp = zeros(self.vector_shape, complex_)
        Ay = zeros_like(y)
        for ii in range(Nv):
            vp += vs[ii]*y[ii]
        Tv = self.equation.rmatvec(vp) - conj(mu)*vp
        for jj in range(Nv):
            Ay[jj] = dot(conj(vs[jj]), Tv)
        return Ay

    def projected_rq(self, vs, mu, x, y):
        ii=0; revchange=inf
        while revchange>1e-12 and ii<self.rayleigh_iterations:
            ii+=1
            self.equation.set_lambda(mu)
            self.jacobian.set_lambda(mu)

            Ax = self.projected_matvec(self.subspace, mu, x)
            Apx = self.projected_matvec(self.subspace, mu, x, prime=True)
            mu_new = mu - dot(conj(y), Ax)/dot(conj(y), Apx)

            muchange = abs(mu_new-mu)
            mu = mu_new
        return mu

    def solve_projected_problem(self, mu, mstart=0, num=1, algorithm=0, xt=None, yt=None, skip_on_unfound=0, tol=1e-6):
        "Solve for the dense NL problem by SLP for m modes"
        mus = zeros(num, complex_)
        res = zeros(num, float_)
        xs = zeros((num,self.Nsubspace), complex_)
        ys = zeros((num,self.Nsubspace), complex_)
        
        #Solve for Nv eigenvalues
        k=0; offset=mstart; mu_next=x_next=y_next=0
        while (k<num) and (offset-mstart<5):
            #Inner SLP loop
            if algorithm==0:
                mu_ip,mu_next,x,y,x_next,y_next,err = self.inner_solver_slp(mu, m=k+offset, tol=tol)
            elif algorithm==1:
                mu_ip,x,y,err = self.inner_solver_ri(mu, xt=xt, yt=yt, tol=tol)
            elif algorithm==2:
                mu_ip,mu_next,x,y,x_next,y_next,err = self.inner_solver_slp(mu, m=0, ordermethod=1, tol=tol)
            elif algorithm==3:
                mu_ip,x,y,err = self.inner_solver_linear(mu, m=k+offset, ordermethod=1)

            #Add converged eigensolution to list
            if err>tol:
                if skip_on_unfound==1:
                    offset+=1
                elif skip_on_unfound==2:
                    mu_ip,x,y,err = self.inner_solver_linear(mu_ip, m=k+offset, ordermethod=1)
                else:
                    err=0

            if err<tol:
                mus[k], xs[k], ys[k], res[k], = mu_ip, x, y, err
                k+=1
        
        if err>tol: offset=-1
        return mus, mu_next, xs, ys, x_next, y_next, res, offset-mstart

    def inner_solver_linear(self, mu, m=0, ordermethod=0):
        A,Ap = self.build_projected_problem(mu)
        thetas, evl, evr = splinalg.eig(A, Ap, left=True, right=True)

        res = 0

        #Iterate on mth ordered theta
        if ordermethod==0:
            order = (real(thetas)+abs(imag(thetas))).argsort()
        elif ordermethod==1:
            order = abs(thetas).argsort()

        mu = mu - thetas[order[m]]
        x = evr[:,order[m]]
        y = evl[:,order[m]]

        logging.info("[LIN] Best match %s" % (mu))
        return mu, x, y, 0

    def inner_solver_slp(self, mu, m=0, ordermethod=0, tol=1e-6):
        "Solve for the dense NL problem by SLP for m modes"
        deltamu = inf; niter = 0
        while (abs(deltamu)>tol) and (niter<20):
            A,Ap = self.build_projected_problem(mu)
            thetas, evl, evr = splinalg.eig(A, Ap, left=True, right=True)

            #Iterate on mth ordered theta
            if ordermethod==1:
                order = abs(thetas).argsort()
            else:
                order = (real(thetas)+abs(imag(thetas))).argsort()

            #Update mu
            deltamu = thetas[order[m]]
            mu = mu - deltamu
            niter += 1

        if niter>=20:
            print "[SLP] didn't converge"
            print m, mu-thetas[order[:5]]

        #Converged eigenvector
        x = evr[:,order[m]]
        y = evl[:,order[m]]

        order = real(thetas).argsort()
        current_num = argmin(abs(thetas[order]))
        next_num = order[min(current_num+1,len(order)-1)]
        mu_next = mu-thetas[next_num]
        x_next = evr[:,next_num]
        y_next = evl[:,next_num]
        
        logging.info("[SLP] res: %.2g, iter: %d, pos: %d, m: %d" % (abs(deltamu), niter, current_num, m))
        return mu, mu_next, x, y, x_next, y_next, abs(deltamu)

    def inner_solver_ri(self, mu, m=0, xt=None, yt=None, tol=1e-6):
        "Solve for the dense NL problem by residual inverse iterations"
        maxres = inf; niter = 0
        
        #Start vector from SLP iteration
        xt_ = zeros(self.Nsubspace, complex_)
        yt_ = zeros(self.Nsubspace, complex_)
        if xt is None:
            xt_ = random(self.Nsubspace).astype(complex_)
        else:
            xt_[:len(xt)] = xt
        if yt is None:
            yt_ = random(self.Nsubspace).astype(complex_)
        else:
            yt_[:len(yt)] = yt
        
        Asig = self.build_projected_problem(mu, prime=False)
        cAsig = conj(Asig.T)

        xt=xt_; yt=yt_          
        while (abs(maxres)>tol) and (niter<30):
            self.equation.set_lambda(mu)
            Tx = self.projected_matvec(self.subspace, mu, xt)
            Thy = self.projected_rmatvec(self.subspace, mu, yt)
            xres = linalg.solve(Asig, Tx)
            yres = linalg.solve(cAsig, Thy)
            maxres = max(abs(xres).max(), abs(yres).max())
        
            xt -= xres; yt -= yres
            self.normalize_vector(xt); self.normalize_vector(yt)

            #Reshift
            if niter>0 and mod(niter,10)==0:
                Asig = self.build_projected_problem(mu, prime=False)
                cAsig = conj(Asig.T)

            #Rayleigh quotient
            mu_new = self.projected_rq(self.subspace, mu, xt, yt)
            mu = mu_new
            niter += 1
        
        logging.info("[RI] %s to %.2g in %d iterations" % (mu, maxres, niter))
        return mu, xt, yt, abs(maxres)*0

    def construct_ritz_vector(self, vs, x):
        v = zeros(self.vector_shape, complex_)
        for ii in range(len(vs)):
            v += x[ii]*vs[ii]
        return v

    def linear_orthogonalize(self, vs, x):
        "Orthogonalize against previous vectors and add to subspace"
        #Modified Gram-Schmidt orthogonalization
        for vk in vs:
                                    x -= dot(conj(vk),x)*vk/dot(conj(vk),vk)
        self.normalize_vector(x)                            #Normalize
        return x

    def linear_biorthogonalize(self, vs, ws, vnew, wnew):
        "Bi-Orthogonalize against previous vectors"

        #Modified Gram-Schmidt orthogonalization
        for kk in range(self.Nsubspace):
            vk, wk = vs[kk], ws[kk]
            alpha = dot(conj(wk),vk)
            
            wnew -= conj(dot(conj(wnew),vk)/alpha)*wk
            vnew -= dot(conj(wk),vnew)/alpha*vk
        
        return self.binormalize_vector(vnew,wnew)

    def normalize_vector(self,x):
        x /= sqrt(dot(conj(x),x))
        return x

    def binormalize_vector(self, v,w):
        #Fix norm of right vector to be one
        v /= linalg.norm(v)
        w /= linalg.norm(w)
        ngamma = dot(conj(w),v)
        w /= conj(ngamma)
        return v, w
    

    def set_sigma(self, sigma):
        self.update_lambda(sigma)
        self.si.set_shift(self.matrix, complex(sigma))

    def local_restart(self):
        pass

class NLArnoldiSolver(KrylovIteration):

    def check_ip(self, mu1, mu2, y, v):
        self.equation.set_lambda(mu1)
        Tv1 = self.equation.matvec(v)

        self.equation.set_lambda(mu2)
        Tv2 = self.equation.matvec(v)

        ip = dot(conj(u), Tv1-Tv2) - (mu1-mu2)*dot(conj(u), v)
        return ip

    def build_projected_problem(self, mu, prime=True):
        Nv = self.Nsubspace
        vs = self.subspace
        
        A = zeros((Nv,Nv), complex_)
        if prime: Ap = zeros((Nv,Nv), complex_)

        #Create projected problem at the projected evalue mu
        self.equation.set_lambda(mu)
        self.jacobian.set_lambda(mu)
        for ii in range(Nv):
            Tv = self.equation.matvec(vs[ii])
            if prime: Tpv = self.jacobian.matvec(vs[ii])
            for jj in range(Nv):
                A[jj,ii] = dot(conj(vs[jj]), Tv) - mu*dot(conj(vs[jj]), vs[ii])
                if prime: Ap[jj,ii] = dot(conj(vs[jj]), Tpv) - dot(conj(vs[jj]), vs[ii])

        if prime:
            return A, Ap
        else:
            return A

    def initialize_subspace(self, vstart=None, N=5, add_converged=0):
        '''
        Inverse iterations for the linear problem near
        '''
        self.subspace = []
        self.Nsubspace = 0

        if add_converged:
            for m in self.modes:
                self.append_subspace(m.right)

        #Starting vector
        if vstart is None: vstart = random(self.vector_shape).astype(self.dtype)
        v = self.normalize_vector(vstart)

        #Linear Arnoldi iterations
        for ii in range(N):
            v = self.si.matvec(v)
            self.append_subspace(v, orthogonalize=True)

    def append_subspace(self, v, orthogonalize=False, normalize=True):

        if orthogonalize:
            self.linear_orthogonalize(self.subspace, v)
        elif normalize:
            self.normalize_vector(v)

        self.subspace.append(v)                 #Add to subspace
        self.Nsubspace += 1                         #Update subspace size

    def calculate(self, number=inf):
        #Create matrix
        self.create()
        self.si = ShiftInvertBlock(overwrite=False)

        self.vector_shape = prod(self.base_shape+(2,))

        #Starting mu guess
        number = int(min(self.totalnumber, number))
        neffrange = self.bracket
        neffapprox = complex(self.bracket[1]*(1.0-1e-5))

        xt = yt = None
        reshift = True
        offset = 0
        self.convergence = []

        ii=0
        isfinished=0
        while not isfinished and ii<number:
            #Choose new neffapprox if in list mode
            if self.nefflist is not None:
                neffapprox = self.nefflist[ii]
            elif self.modelist is not None:
                neffapprox = self.modelist[ii].neff
            else:
                pass

            mu = neffapprox**2*self.k0**2
            
            if 1:
                #Restart subspace with linear Arnoldi iterations
                sigma = mu
                self.set_sigma(sigma)
                self.initialize_subspace(N=5)

            niter=0; resmax = inf; tol = 1e-10
            while (resmax>tol):
                #Solve projected problem
                ptol = tol if niter>10 else min(resmax, 1e-6)
                mus, mu_next, xs, ys, x_next, y_next, pres, skip = self.solve_projected_problem(mu, ii+offset, algorithm=0, tol=ptol)

                #If projected problem has had difficulty being solved use the RQ
                if skip<0:
                    mu = self.rayleigh_quotient(mu, v)
                    print "[RQ]:", mu
                else:
                    v = self.construct_ritz_vector(self.subspace, xs[-1])
                    evchange = abs(mu-mus[-1]); mu = mus[-1]

                #Residual vector
                self.equation.set_lambda(mu)
                vres = self.equation.matvec(v)-mu*v
                
                resmax = abs(vres).max()

                #Skip over difficult modes if we're sure
                if skip>0 and resmax<1e-3:
                    print "Near convergence, skipping modes"
                    offset += skip

                if resmax>tol:
                    #Residual inverse iteration
                    v = self.si.matvec(vres)

                    #Orthogonalize and add to subspace
                    self.append_subspace(v, orthogonalize=True)

                niter+=1
                self.convergence.append(resmax)
                print " [%d] Subspace: %d, Offset: %d, Max residual: %.2g, EV change; %.2g" \
                    % (ii, self.Nsubspace, offset, resmax, evchange)
            
            #Estimate convergence rate
            xc = array(self.convergence[-min(5,niter):])
            convergence_rate = -mean(diff(log10(xc)))
            print "Converged with final convergence rate:", convergence_rate, "Skip:", skip
            print

            #Store mode
            nextmode = self.mode_class(coord=self.equation.coord, symmetry=self.wg.symmetry, \
                m0=self.m0, wl=self.wl, evalue=mu, wg=self.wg)
            nextmode.residue = (0,resmax)
            nextmode.right = self.normalize_vector(v)
            self.modes.append(nextmode)     

            reshift = False
            
            mu_next = mu
            if ii<(number-1):
                #Choose next eigensolution guess
                mus, xs, ys, pres, skip = self.solve_projected_problem(mu, ii+offset+1, algorithm=3, tol=1e-3)
                v = self.construct_ritz_vector(self.subspace, xs[-1])
                mu_next = mus[-1]; xt = xs[-1]; yt = ys[-1]
                print "Selected next ev:", mu_next, "  Skip:", skip

            #Restart if needed
            if self.Nsubspace>40:
                upper_limit = real(mu)
                sigma = real(mu)*0.999      #Choose new shift
                reshift = True                              #Reshift
                print "Restarted at", sigma

            #Reshift if needed
            if reshift:
                self.set_sigma(sigma)
                self.initialize_subspace(vstart=v, N=5, add_converged=1)
            
            ii += 1
            mu = mu_next


class NLBiArnoldiSolver(KrylovIteration):

    def initialize_subspace(self, mu, vstart=None, wstart=None, N=5):
        '''
        Inverse iterations for the linear problem near
        '''
        self.lsubspace = []
        self.rsubspace = []
        self.Nsubspace = 0

        #starting vector
        if vstart is None: vstart = random(self.vector_shape).astype(self.dtype)
        if wstart is None: wstart = random(self.vector_shape).astype(self.dtype)

        v,w = self.binormalize_vector(vstart,wstart)

        for ii in range(N):
            v = self.si.matvec(v)                                           #Inverse iteration
            w = self.si.rmatvec(w)                                          #Inverse iteration
            self.append_subspace(v,w, mu, orthogonalize=True)

    def append_subspace(self, v, w, mu, orthogonalize=False, normalize=True):
        if orthogonalize:
            self.linear_biorthogonalize(self.rsubspace,self.lsubspace,v,w)

        if self.Northogonalize>0:
            self.approximate_orthogonalize(v,w, mu)

        if normalize:
            self.binormalize_vector(v,w)

        self.rsubspace.append(v)                #Add to subspace
        self.lsubspace.append(w)
        self.Nsubspace += 1                     #Update subspace size

        self.augment_projected_problem(self.Nsubspace-1)
        
    @timer.time_function()
    def build_projected_problem_old(self, mu, prime=True):
        Nv = self.Nsubspace
        A = zeros((Nv,Nv), complex_)
        if prime: Ap = zeros((Nv,Nv), complex_)

        vs = self.rsubspace; ws = self.lsubspace

        #Create projected problem at the projected evalue mu
        self.equation.set_lambda(mu)
        self.jacobian.set_lambda(mu)
        for ii in range(Nv):
            Tv = self.equation.matvec(vs[ii])
            if prime: Tpv = self.jacobian.matvec(vs[ii])
            for jj in range(Nv):
                A[jj,ii] = dot(conj(ws[jj]), Tv) - mu*dot(conj(ws[jj]), vs[ii])
                if prime: Ap[jj,ii] = dot(conj(ws[jj]), Tpv) - dot(conj(ws[jj]), vs[ii])

        if prime:
            return A, Ap
        else:
            return A

    @timer.time_function()
    def augment_projected_problem(self, ii):
        #Augment the constant part of the projected problem
        #with the new vectors:
        Nv = self.Nsubspace
        self.equation.set_lambda(0, onlyconst=1)
        vs = self.rsubspace; ws = self.lsubspace

        Tv = self.equation.matvec(vs[ii])
        Tu = self.equation.rmatvec(ws[ii])
        for jj in range(Nv):
            self.pp_const[jj,ii] = dot(conj(ws[jj]), Tv)
            self.pp_linear[jj,ii] = dot(conj(ws[jj]), vs[ii])
            
            self.pp_const[ii,jj] = dot(conj(Tu), vs[jj])
            self.pp_linear[ii,jj] = dot(conj(ws[ii]), vs[jj])

    @timer.time_function()
    def build_projected_problem_new(self, mu, prime=True):
        Nv = self.Nsubspace
        A = self.pp_const[:Nv,:Nv] - mu*self.pp_linear[:Nv,:Nv]
        if prime: Ap = -self.pp_linear[:Nv,:Nv]+0

        vs = self.rsubspace; ws = self.lsubspace

        #Create projected problem at the projected evalue mu
        self.equationip.set_lambda(mu, onlyev=1)
        self.jacobian.set_lambda(mu)
        for ii in range(Nv):
            Tv = self.equationip.matvec(vs[ii])
            if prime: Tpv = self.jacobian.matvec(vs[ii])
            for jj in range(Nv):
                A[jj,ii] += dot(conj(ws[jj]), Tv)
                if prime: Ap[jj,ii] += dot(conj(ws[jj]), Tpv)

        if prime:
            return A, Ap
        else:
            return A

    build_projected_problem=build_projected_problem_new

    @timer.time_function()
    def calculate(self, number=inf):
        self.create()
        self.si = ShiftInvertBlock(overwrite=False)
        self.vector_shape = prod(self.base_shape+(2,))

        Nvmax = 80
        self.pp_const = zeros((Nvmax,Nvmax), complex_)
        self.pp_linear = zeros((Nvmax,Nvmax), complex_)
        self.rorthogonalize = []
        self.lorthogonalize = []
        self.Northogonalize = 0

        #Starting mu guess
        number = int(min(self.totalnumber, number))
        neffrange = self.bracket
        neffapprox = complex(self.bracket[1]*(1.0-1e-5))
        mu_next = neffapprox**2*self.k0**2
        v_next = random(self.vector_shape).astype(self.dtype)
        w_next = random(self.vector_shape).astype(self.dtype)

        xt = yt = None
        reshift = True
        offset = 0
        self.convergence = []

        ii=0
        isfinished=0
        while not isfinished and ii<number:
            #Choose new neffapprox if in list mode
            if self.nefflist is not None:
                neffapprox = self.nefflist[ii]
                mu = neffapprox**2*self.k0**2
                v = random(self.vector_shape).astype(self.dtype)
                w = random(self.vector_shape).astype(self.dtype)
                
                #vx = v.reshape(self.base_shape+(2,))
                #vx *= exp(-self.equation.coord.ms[:,newaxis]**2)
                #wx = w.reshape(self.base_shape+(2,))
                #wx *= exp(-self.equation.coord.ms[:,newaxis]**2)

            elif self.modelist is not None:
                neffapprox = self.modelist[ii].neff
                mu = neffapprox**2*self.k0**2
                v = self.modelist[ii].right + 0
                w = self.modelist[ii].left + 0

            else:
                mu = mu_next
                w = w_next
                v = v_next

            restart = True
            if mod(ii,40)==0:
                restart = True

            if restart:
                print "Restart at neff=", sqrt(mu)/self.k0
                #Restart subspace with linear Arnoldi iterations
                sigma = mu
                mu = mu
                self.set_sigma(sigma)
                self.initialize_subspace(mu, wstart=w, vstart=v, N=5)

            niter=0; resmax = inf; tol = 1e-10
            while (resmax>tol):
                #Solve projected problem
                ptol = tol if niter>10 else min(resmax, 1e-6)
                mus, mu_next, xs, ys, x_next, y_next, pres, skip = self.solve_projected_problem(mu, ii+offset, algorithm=2, tol=ptol)

                #If projected problem has had difficulty being solved use the RQ
                if skip<0:
                    mu = self.rayleigh_quotient(mu, v)
                    print "[RQ]:", mu
                else:
                    v = self.construct_ritz_vector(self.rsubspace, xs[-1])
                    w = self.construct_ritz_vector(self.lsubspace, ys[-1])
                    evchange = abs(mu-mus[-1]); mu = mus[-1]

                #Residual vector
                self.equation.set_lambda(mu)
                vres = self.equation.matvec(v)-mu*v
                wres = self.equation.rmatvec(w)-conj(mu)*w
        
                resmax = max(abs(vres).max(), abs(wres).max())
                if resmax>tol:
                    #Residual inverse iteration
                    v = self.si.matvec(vres)
                    w = self.si.rmatvec(wres)

                    #Orthogonalize and add to subspace
                    self.append_subspace(v,w,mu, orthogonalize=True)

                niter+=1
                self.convergence.append(resmax)
                print "[%d] Approximate eigenvalue: %s" % (niter, mu)
                print " Subspace: %d, Offset: %d, Max residual: %.2g, EV change; %.2g" \
                    % (self.Nsubspace, offset, resmax, evchange)
                print
    
            #Estimate convergence rate
            xc = array(self.convergence[-min(5,niter):])
            convergence_rate = -mean(diff(log10(xc)))
            print "Converged with final convergence rate:", convergence_rate, "Skip:", skip
            print "Closest NL evalue:", mu_next
            
            #mus,_,xs,ys,_,_,_,_ = self.solve_projected_problem(mu, 0, num=5, algorithm=3, tol=ptol)
            #print "Closest L evalues:", mus

            #Store vector
            nextmode = self.mode_class(coord=self.equation.coord, symmetry=self.wg.symmetry, \
                m0=self.m0, wl=self.wl, evalue=mu, wg=self.wg)
            nextmode.residue = (0,resmax)
            nextmode.right = self.normalize_vector(v)
            nextmode.left = self.normalize_vector(w)
            self.modes.append(nextmode)

            #Next guess
            v_next = self.construct_ritz_vector(self.rsubspace, x_next)
            w_next = self.construct_ritz_vector(self.lsubspace, y_next)

            #Add converged vector to orthogonalized list
            self.rorthogonalize.append(self.normalize_vector(v))
            self.lorthogonalize.append(self.normalize_vector(w))
            self.Northogonalize += 1
            ii+=1
        return self.modes       

    def check_ip(self, mu1, mu2, y, v):
        self.equation.set_lambda(mu1)
        Tv1 = self.equation.matvec(v)

        self.equation.set_lambda(mu2)
        Tv2 = self.equation.matvec(v)

        ip = dot(conj(u), Tv1-Tv2) - (mu1-mu2)*dot(conj(u), v)
        return ip


DefaultSolver = NLArnoldiSolver

