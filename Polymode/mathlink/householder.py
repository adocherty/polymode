# coding=utf-8
'''
┈┈┈┈┈┈┈┈┈┈ Householder Transformations  ┈┈┈┈┈┈┈┈┈┈
Elementary Householder transformations
QR decomposition implemented with Householder transformations

┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈┈
'''
from __future__ import division

from numpy import *
import os, sys, logging

#
# Helpful elementary transforms
#
class ElementaryTransform(object):
    def apply(self, X):
        pass
    def rapply(self, X):
        pass
    def inverse(self, X):
        pass
    def rinverse(self, X):
        pass
    def matvec(self, x):
        Qx = x.copy()
        self.apply(Qx)
        return Qx

class HouseholderTransform(ElementaryTransform):
    "Defines and applied implicit householder transform"
    def __init__(self, v, offset='end'):
        #Ensure that v is not divided by zero
        vnorm = sqrt(dot(conj(v),v))
        if vnorm>0:
            self.v = v/vnorm
        else:
            self.v = v
        self.set_range(offset)
    
    def __str__(self):
        return "Householder length %d" % (len(self.v),)
    def __repr__(self):
        return "Householder length %d" % (len(self.v),)
    
    def set_range(self, offset):
        """
        Default action is last N elements
        use offset to change this behaviour
        """
        if offset.lower() == 'start':
            self.range = slice(None, len(self.v))
        elif offset.lower() == 'end':
            self.range = slice(-len(self.v), None)
    
    def apply(self, X):
        """
        Apply householder transformation to columns of X
        so Xnew = Q X
        """
        #Vector apply:
        #for x in atleast_2d(X.T):
        #   x[self.range] -= 2*dot(conj(self.v), x[self.range])*self.v
        Y = atleast_2d(X)[self.range]
        Y -= 2*dot(self.v[:,newaxis], dot(conj(self.v),Y)[newaxis,:])


    #Inverse is the same as action
    inverse = apply

    def rapply(self, X):
        """
        Apply right householder transformation to rows of X
        so Xnew = X Q
        """
        #For vector:
        #for x in atleast_2d(X):
        #   x[self.range] -= 2*dot(x[self.range], self.v)*conj(self.v)
        Y = atleast_2d(X)[:,self.range]
        Y -= 2*dot(dot(Y, self.v)[:,newaxis], conj(self.v)[newaxis,:])
        
    #Inverse is the same as action
    rinverse = rapply

class UnitHouseholderTransform(HouseholderTransform):
    def __init__(self, v, k=0, offset='end'):
        #Take alpha = - exp(i arg(x_k)) ||v||
        #Note exp(i arg(v_k)) = v_k/abs(v_k) only if v_k<>0
        ve1 = v.copy()
        self.alpha = -sqrt(dot(conj(v),v)) * exp(1j*angle(v[k]))
        ve1[k] -= self.alpha

        HouseholderTransform.__init__(self, ve1, offset)

class GivensRotation(ElementaryTransform):
    "Defines and applied implicit householder transform"
    def __init__(self, cs=(0,0), locations=(0,0)):
        self.c, self.s = cs
        self.locations = locations
    
    def set_offset(self, offset):
        """
        Default action is to count i,j from the start
        use offset to change this behaviour
        """
        self.offset

    def rotate(self, X, a, b):
        #**todo** This is damn slow .. consider a speed up!
        i,j = self.locations
        for x in atleast_2d(X):
            xn = a*x[i] + b*x[j]
            yn = -b*x[i] + a*x[j]
            x[i]=xn; x[j] = yn
    
    def apply(self, X):
        """
        Apply Givens rotation to columns of X
        so Xnew = G X
        """
        self.rotate(X.T, self.c, self.s)

    def inverse(self, X):
        """
        Apply inverse Givens rotation to columns of X
        so Xnew = G^-1 X
        """
        self.rotate(X.T, self.c, -self.s)

    def rapply(self, X):
        """
        Apply right Givens rotation to rows of X
        so Xnew = X G
        """
        self.rotate(X, self.c, -self.s)

    def rinverse(self, X):
        """
        Apply inverse right Givens rotation to rows of X
        so Xnew = X G
        """
        self.rotate(X, self.c, self.s)
        
class ZeroGivensRotation(GivensRotation):
    """
    A Givens rotation defined to introduce a zero to v[j]
    """
    def __init__(self, v, locations=(0,0)):
        i,j = locations
        
        #Handle different cases of coefficients.
        #See:
        # http://en.wikipedia.org/wiki/Givens_rotation#Stable_calculation
        #
        # and
        #
        #Anderson, Edward (2000)
        #Discontinuous Plane Rotations and the Symmetric Eigenvalue Problem
        #LAPACK Working Note 150, University of Tennessee, UT-CS-00-454, December 4, 2000.
        #http://www.netlib.org/lapack/lawnspdf/lawn150.pdf
        #
        #D. Bindel, J. Demmel, W. Kahan, O. Marques. (2001)
        #On Computing Givens rotations reliably and efficiently
        #LAPACK Working Note 148, University of Tennessee, UT-CS-00-449, January 31, 2001.
        #http://www.netlib.org/lapack/lawnspdf/lawn148.pdf
        #
        if v[i]==0 and v[j]==0:
            c = 1; s = 0; r = 0
            
        elif v[j]==0:
            c = sign(v[i]); s = 0; r = abs(v[i])

        elif v[i]==0:
            c = 0; s = sign(v[j]); r = abs(v[j])

        elif abs(v[j])>abs(v[i]):
            t = v[i]/v[j]
            u = sign(v[j]) * sqrt(1+t*t)
            s = 1/u
            c = s*t
            r = v[j]*u

        else:
            t = v[j]/v[i]
            u = sign(v[i]) * sqrt(1+t*t)
            c = 1/u
            s = c*t
            r = v[i]*u

        print c,s

        self.r = r
        GivensRotation.__init__(self, cs=(c,s), locations=locations)

class CompoundTransform(object):
    def __init__(self, N=0):
        self.transforms = [None]*N

    def clear_transforms(self):
        self.transforms = []
    
    def set_transform(self, trans, j):
        if len(self.transforms)<=j:
            self.transforms.append( [None]*(j+1-len(self.transforms)))
        self.transforms[j] = trans

    def add_transform(self, trans, offset=None):
        if offset is not None:
            trans.set_range(offset)
        self.transforms.append(trans)

    def add_transforms(self, trans, offset=None):
        for t in trans:
            self.add_transform(t, offset)
    
    def q_action(self, x):
        for trans in self.transforms[::-1]:
            if trans is not None:
                trans.apply( x )

    def q_raction(self, x):
        for trans in self.transforms:
            if trans is not None:
                trans.rapply( x )

    def q_inverse(self, x):
        for trans in self.transforms:
            if trans is not None:
                trans.inverse( x )

    def q_rinverse(self, x):
        for trans in self.transforms[::-1]:
            if trans is not None:
                trans.rinverse( x )

class HouseholderQR(CompoundTransform):
    '''
    Computes Householder QR decompostion and returns R
    '''
    def __init__(self, overwrite = False, offset=0):
        self.overwrite = overwrite
        self.offset = offset
        CompoundTransform.__init__(self)
    
    def __call__(self, V, R=None, jstart=0, jend=None):
        N = V.shape[0]
        if jend is None: jend = N
        
        if self.overwrite:
            Vtmp = V[jstart:jend]
        else:
            Vtmp = V[jstart:jend].copy()
            #Apply prior transforms to all vectors
            #self.q_inverse( Vtmp.T )

        for jj in range(jstart, jend):
            vi = Vtmp[jj-jstart,jj:]

            #Create householder transformation and add to internal list
            hhtrans = UnitHouseholderTransform(vi)
            self.add_transform(hhtrans)

            #Update R matrix, with offset (to allow for Arnoldi updates)
            if jj>=self.offset and R is not None:
                R[:jj,jj-self.offset] = Vtmp[jj-jstart,:jj]
                R[jj,jj-self.offset] = hhtrans.alpha
            
            ## Apply householder reflection to all subsequent vectors
            hhtrans.apply( Vtmp[jj-jstart+1:jend].T )
        return R

class HouseholderHessenberg(CompoundTransform):
    '''
    Computes Householder Hessenburg decomposition
     A = Q* H Q
    '''
    def __init__(self, overwrite = False):
        self.overwrite = overwrite
        CompoundTransform.__init__(self)

    def __call__(self, A):
        N = A.shape[0]

        if self.overwrite: H = A
        else: H = A.copy()

        for jj in range(1, N):
            #Create householder transformation and add to internal list
            hhtrans = UnitHouseholderTransform( H[jj:,jj-1], offset='end' )
            self.add_transform(hhtrans)

            #Update H matrix
            H[jj,jj-1] = hhtrans.alpha
            H[jj+1:,jj-1] = 0
            
            ## Apply householder transformation to partitioned matrix
            hhtrans.apply(H[jj:,jj:])
            hhtrans.rapply(H[:,jj:])

        return H

class RowHessenberg(CompoundTransform):
    '''
    Computes Householder Hessenburg decomposition
     A = Q* H Q
    Along the rows of the matrix, so A[-1,-1] is undisturbed
    '''
    def __init__(self, overwrite = False):
        self.overwrite = overwrite
        CompoundTransform.__init__(self)

    def __call__(self, A):
        N = A.shape[0]

        if self.overwrite: H = A
        else: H = A.copy()

        for jj in range(1, N):
            #Create householder transformation and add to internal list
            hhtrans = UnitHouseholderTransform( H[-jj, :-jj], k=-1, offset='start' )
            self.add_transform(hhtrans)
            
            #Update -jj th row of H matrix
            H[-jj,:-jj] = 0
            H[-jj,-jj-1] = hhtrans.alpha
            
            ## Apply householder transformation to _rows_ of partitioned matrix
            hhtrans.apply(H[:-jj].T)
            hhtrans.rapply(H.T)
        return H

class HessenbergReduce(CompoundTransform):
    """
    Given A in hessenberg form reduce to triangular form
    """
    def __init__(self, overwrite = False):
        self.overwrite = overwrite
        CompoundTransform.__init__(self)

    def __call__(self, A):
        Na = A.shape[0]
        for ii in range(Na-1):
            #Givens transforms for all off-diagonal elements
            givens = ZeroGivensRotation(A[:,ii], (ii,ii+1))
            givens.apply(A)
            
            #To replicate action add in reverse order
            self.transforms = [givens] + self.transforms
        return A
        

