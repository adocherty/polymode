#Cauchy root finding
from __future__ import division

import numpy as np
import pylab as pl
import numpy.linalg as la

from numpy import pi, exp, arange, zeros, ones, real, imag, dot, roots, absolute, inf, roll, log, unwrap, linalg
from numpy import array, angle

_CROOT_PLOTRESULTS=0
_CROOT_DEBUG=0

def findzero_delves(f, fprime, z0=0, R=1, N=None, alpha=1, trange=None, tol=1e-6):
    '''
    Cauchy integral method for finding the zeros of an analytic function
    f: function of a single variable returns function value
    fprime: derivative of function f
    z0: center location in the complex plane
    R: radius of region in which to bound the search
    N: Number of boundary integral points, otherwise automatic

    Algorithim from Delves and Lyness
    '''

    Nt = 128 if N is None else N
    trange = (0,2*pi) if trange is None else trange

    residue = 1.0; niter = 0; maxiter = 5
    while (residue>tol) and (niter<maxiter):
        #Evaluate function on circle radius R
        dt = 2*pi/Nt
        thetas = np.arange(trange[0], trange[1], dt)

        #Calculate Phi'(t)/Phi(t)
        zs = R*np.exp(1j*thetas)
        phi = fprime(zs+z0)/f(zs+z0)*1j*zs/alpha

        if _CROOT_PLOTRESULTS>1:
            pl.plot(thetas, np.unwrap(np.angle(phi)), 'r--')

        #Estimate number of zeros
        I0 = np.real(1/(2*pi*1j)*np.sum(phi)*dt)
        K = int(round(I0))

        #Reject K too large or too small
        if (I0<0.999) or (I0>50):
            print "Warning, no roots found", I0
            return array([])

        #Evaluate integral by trapezoidal rule
        I = np.zeros(K, np.complex_)
        for k in range(K):
            I[k] = 1/(2*pi*1j)*np.sum(zs**(k+1)*(phi))*dt

        #Solve for the coefficients
        ac = np.zeros(K+1, np.complex_)
        ac[0] = 1.0
        for k in range(K):
            ac[k+1] = - np.dot(ac[k::-1], I[:k+1])/(k+1)

        calc_roots = np.roots(ac)

        if _CROOT_PLOTRESULTS>0:
            pl.plot(calc_roots.real, calc_roots.imag, 'b.')

        #Check error
        residue = np.absolute(f(calc_roots+z0)).max()

        #Increase resolution
        Nt = 2*Nt
        niter += 1

    print "Calculated %d roots in %d iterations to a residue %.2g" % (K, niter, residue)
    return calc_roots+z0


def findzero_adr(f, z0=0, R=1, N=None, nroot=None, tol=1e-6, alpha=1, quiet=False, trange=None, maxiter=10):
    '''
    Cauchy integral method for finding the zeros of an analytic function
    doesn't require the derivative of the function
    f: function of a single variable returns function value
    z0: center location in the complex plane
    R: radius of region in which to bound the search
    N: Number of boundary integral points, otherwise automatic

    Algorithm from:
    "A numerical method for locating the zeros and poles of a meromorphic function"
    LF Abd-Ellal, LM Delves, JK Reid - Numerical Methods for Nonlinear Algebraic Equations, 1970
    '''
    trange = (0,2*pi) if trange is None else trange
    Nt = 32 if N is None else N
    T = trange[1]-trange[0]

    residue = inf; niter = 0
    while (residue>tol) and (niter<maxiter):
        #Evaluate function on circle radius R
        dt = T/Nt
        thetas = np.arange(trange[0],trange[1],dt)

        #Careful to 'unwrap' the phase
        zs = R*np.exp(1j*thetas)
        fz = f(zs+z0)
        #fz = abs(fz)*exp(1j*unwrap(angle(fz)/alpha)*alpha)

        if _CROOT_PLOTRESULTS>1:
            pl.plot(thetas, np.unwrap(np.angle(1/fz)), 'r--')

        #Check for zeros close to the boundary

        #Number of roots enclosed
        if nroot is None:
            I0 = (np.unwrap(np.angle(fz))[-1]-np.unwrap(np.angle(fz))[0])/(2*pi)
            K = int(round(I0))
        else:
            I0 = K = nroot

        #Reject K too large or too small
        if (I0<0.99) or (I0>50):
            if not quiet: print "Warning, no roots found", I0
            return array([])

        #Construct companion matrix for the polynomial equation
        #and the truncated matrix of Newton's equations.
        Ic = np.zeros(K, complex)
        XC = np.zeros((K,K), complex)
        for k in range(K):
            Ic[k] = R**(k+K+1)*np.sum(np.exp(1j*(k+K+1)*thetas)/fz)*dt
            for m in range(K):
                XC[k,m] = R**(k+m+1)*np.sum(np.exp(1j*(k+m+1)*thetas)/fz)*dt

        #Solve for the coefficients
        ac = np.ones(K+1, complex)
        ac[1:] = la.solve(XC,-Ic)[::-1]

        calc_roots = np.roots(ac)

        if _CROOT_PLOTRESULTS>0:
            pl.plot(calc_roots.real, calc_roots.imag, 'kx')

        #Check error
        residue = absolute(f(calc_roots+z0)).max()

        #Increase resolution
        Nt = 2*Nt
        niter += 1

    if not quiet: print "Calculated %d roots in %d iterations with approximate error %.2g" % (K, niter, residue)
    return calc_roots+z0

def findzero_carpentier(f, z0=0, R=1, N=None, tol=1e-6, alpha=1, trange=None, quiet=False, force=False, maxiter=10):
    '''
    Cauchy integral method for finding the zeros of an analytic function
    doesn't require the derivative of the function
    f: function of a single variable returns function value
    z0: center location in the complex plane
    R: radius of region in which to bound the search
    N: Number of boundary integral points, otherwise automatic

    Algorithim from Carpentier and dos Santos
    '''
    trange = (0,2*pi) if trange is None else trange
    Nt = 32 if N is None else N
    T = trange[1]-trange[0]

    residue = inf; niter = 0
    while (residue>tol) and (niter<maxiter):

        #Evaluate function on circle radius R
        dt = T/Nt
        thetas = np.arange(trange[0],trange[1],dt)

        zs = R*np.exp(1j*thetas)
        zshift = R*np.exp(1j*(thetas-dt))

        #Careful to 'unwrap' the phase of the root
        fz = f(zs+z0)
        #fz_dt = f(zshift+z0)
        fz_dt = np.roll(fz,1)

        #Take the correct branch of the log function
        g = np.log(fz/fz_dt)
        g = np.real(g) + np.unwrap(np.imag(g))*1j
        
        if _CROOT_PLOTRESULTS==2:
            pl.plot(np.real(zs+z0)/f.k0, np.abs(fz), 'b--')

        if _CROOT_PLOTRESULTS==3:
            pl.plot(thetas, np.unwrap(np.angle(g)), 'b--')
            pl.plot(thetas, np.angle(fz), 'g-')

        #Check for zeros close to the boundary

        #Number of roots enclosed
        I0 = np.real(np.sum(g)/(T*1j))/alpha

        if not quiet: print "Roots found:", I0

        if not np.isfinite(I0) or (not force and (I0<0.999 or I0>50)):
            if not quiet: print "No roots were found."
            return array([])

        K = int(round(I0))
        
        #Calculate the contour integrals
        Ic = np.zeros(K+1, complex)
        for k in range(1,K+1):
            Ic[k] = (R**k)*np.sum(g*np.exp(1j*k*thetas))*k \
                     /(np.exp(1j*k*T/Nt)-1)/Nt/alpha

        #Construct companion matrix for the polynomial equation
        #and the truncated matrix of Newton's equations.
        XC = np.zeros((K,K), np.complex_)
        X = np.zeros((K,K), np.complex_)
        for k in range(0,K):
            X[k,k] = K-k
            XC[k,k:] = Ic[1:K+1-k]
            if k>0:
                X[k-1,k:] = Ic[1:K-k+1]
                XC[k,k-1] = K-k

        #Find eigenvalues - the roots of the equation
        calc_roots = la.eigvals(dot(XC,linalg.inv(X)))
        #calc_roots = linalg.eigvals(dot(linalg.inv(X),XC))

        if _CROOT_PLOTRESULTS>0:
            pl.plot(calc_roots.real, calc_roots.imag, 'kx')

        #Check error
        residue = np.absolute(f(calc_roots+z0)).max()

        if not quiet: print "Roots:", calc_roots+z0
        if not quiet: print "Res:", residue, "at N=", Nt
        
        #Increase resolution
        Nt = 2*Nt
        niter += 1

    if not quiet:
        print "Calculated %d roots in %d iterations with approximate error %.2g" % (K, niter, residue)
    return calc_roots+z0


def findzero_matrix(A, z0=0, R=1, N=None, tol=1e-6, maxiter=10):
    '''
    Cauchy integral method for finding the zeros of an analytic matrix function
    doesn't require the derivative of the function
    A: matrix function of a single variable
    z0: center location in the complex plane
    R: radius of region in which to bound the search
    N: Number of boundary integral points, otherwise automatic

    Algorithim from "Foundations of Photonic Crystal Fibers", Zolla et al
    '''

    Nt = 64 if N is None else N
    residue = 1.0
    niter = 0
    maxiter = 5
    while (residue>tol) and (niter<maxiter):

        #Evaluate function on circle radius R
        dt = 2*pi/Nt
        thetas = arange(0,2*pi,dt)

        zs = R*exp(1j*thetas)
        fz = zeros(A(0).shape+zs.shape, complex_)
        for ii in range(Nt):
            try:
                fz[...,ii] = inv(A(zs[ii]+z0))
            except:
                break

        #Matrix Cauchy integrals
        I1 = 1/(2*pi*1j)*sum(fz*1j*zs, axis=-1)*dt
        I2 = 1/(2*pi*1j)*sum(zs*fz*1j*zs, axis=-1)*dt

        #Diagonalize
        v,w = linalg.eig(I1)

        calc_roots = dot(linalg.inv(w),dot(I2,w)).diagonal()/v
        residue = absolute([linalg.det(A(z+z0)) for z in calc_roots]).max()

        #Increase resolution
        Nt = 2*Nt
        niter += 1

    print "Calculated x roots in %d iterations with approximate error %.2g" % (niter, residue)
    return calc_roots+z0
