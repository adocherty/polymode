# Cython interface for amos

import numpy as np
cimport numpy as np

#The python complex double
cdef extern from "complexobject.h":
    ctypedef struct Py_complex:
        double real
        double imag

    ctypedef class __builtin__.complex [object PyComplexObject]:
        cdef Py_complex cval

#Wrap the complex double class
cdef extern from "complex":
    ctypedef struct c_double "std::complex<double>":
        double real()
        double imag()
    c_double *new_ccomplex "new std::complex<double>" (double x, double y)
    void del_complex "delete" (c_double *c)

ctypedef np.int_t ITYPE_t
ctypedef np.float_t DTYPE_t

#My external defs
cdef extern from "amos.h":
    c_double c_besselj "amos::besselJ"(double nu, c_double x)
    c_double c_hankel "amos::hankel"(int m, double nu, c_double x)
    c_double c_hankelp "amos::hankelp"(int m, double nu, c_double x)

cdef extern from "bessel_function_ratios.h":
    c_double c_bessel_ratio "besselratio::besselj_ratio"(double nu, c_double x)
    c_double c_hankel_ratio "besselratio::hankel1_ratio"(double nu, c_double x)

def besselj(double m, complex x):
    cdef c_double *cxp = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cy = c_besselj(m, cxp[0])
    return complex(cy.real(),cy.imag())

#The outward hankel function
def hankel1(double m, complex x):
    cdef c_double *cxp = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cy = c_hankel(1, m, cxp[0])
    return complex(cy.real(),cy.imag())

def hankel1p(double m, complex x):
    cdef c_double *cxp = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cy = c_hankelp(1, m, cxp[0])
    return complex(cy.real(),cy.imag())

#The inward hankel function
def hankel2(double m, complex x):
    cdef c_double *cxp = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cy = c_hankel(2, m, cxp[0])
    return complex(cy.real(),cy.imag())

def hankel2p(double m, complex x):
    cdef c_double *cxp = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cy = c_hankelp(2, m, cxp[0])
    return complex(cy.real(),cy.imag())


#The Bessel ratio function J_{m+1}(x)/J_{m-1}(x)
#Currently:
# * m MUST be float64
# * ans MUST be complex128

def _besselj_ratio_1d(np.ndarray[np.float_t, ndim=1] m, complex x, np.ndarray[np.cdouble_t, ndim=1] ans):
    cdef int N1 = m.shape[0]
    cdef c_double *cx_p = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cans

    cdef int jj, kk
    for jj in range(N1):
        cans = c_bessel_ratio(m[jj], cx_p[0])
        ans[jj].real = cans.real()
        ans[jj].imag = cans.imag()

def _besselj_ratio_2d(np.ndarray[np.float_t, ndim=2] m, complex x, np.ndarray[np.cdouble_t, ndim=2] ans):
    cdef int N1 = m.shape[0], N2 = m.shape[1]
    cdef c_double *cx_p = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cans

    cdef int jj, kk
    for jj in range(N1):
        for kk in range(N2):
            cans = c_bessel_ratio(m[jj,kk], cx_p[0])
            ans[jj,kk].real = cans.real()
            ans[jj,kk].imag = cans.imag()

def _besselj_ratio(double m, complex x):
    cdef c_double *cxp = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cy = c_bessel_ratio(m, cxp[0])
    return complex(cy.real(),cy.imag())

def besselj_ratio(m, x):
    "Calculates the ratio of Bessel functions"
    ma = np.array(m).astype(np.float64)
    x = complex(x)

    if ma.ndim==0:
        br = _besselj_ratio(ma, x)
        
    elif ma.ndim==1:
        br = np.zeros(m.shape, dtype=np.complex128)
        _besselj_ratio_1d(ma, x, br)

    elif ma.ndim==2:
        br = np.zeros(m.shape, dtype=np.complex128)
        _besselj_ratio_2d(ma, x, br)

    return br

#The Bessel ratio function J_{m+1}(x)/J_{m-1}(x)
#Currently:
# * m MUST be float64
# * ans MUST be complex128

def _hankel1_ratio_1d(np.ndarray[np.float_t, ndim=1] m, complex x, np.ndarray[np.cdouble_t, ndim=1] ans):
    cdef int N1 = m.shape[0]
    cdef c_double *cx_p = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cans

    cdef int jj, kk
    for jj in range(N1):
        cans = c_hankel_ratio(m[jj], cx_p[0])
        ans[jj].real = cans.real()
        ans[jj].imag = cans.imag()

def _hankel1_ratio_2d(np.ndarray[np.float_t, ndim=2] m, complex x, np.ndarray[np.cdouble_t, ndim=2] ans):
    cdef int N1 = m.shape[0], N2 = m.shape[1]
    cdef c_double *cx_p = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cans

    cdef int jj, kk
    for jj in range(N1):
        for kk in range(N2):
            cans = c_hankel_ratio(m[jj,kk], cx_p[0])
            ans[jj,kk].real = cans.real()
            ans[jj,kk].imag = cans.imag()

def _hankel1_ratio(double m, complex x):
    cdef c_double *cxp = new_ccomplex(x.cval.real, x.cval.imag)
    cdef c_double cy = c_hankel_ratio(m, cxp[0])
    return complex(cy.real(),cy.imag())

def hankel1_ratio(m, x):
    "Calculates the ratio of Bessel functions"
    ma = np.array(m).astype(np.float64)
    x = complex(x)

    if ma.ndim==0:
        br = _hankel1_ratio(ma, x)
        
    elif ma.ndim==1:
        br = np.zeros(m.shape, dtype=np.complex128)
        _hankel1_ratio_1d(ma, x, br)

    elif ma.ndim==2:
        br = np.zeros(m.shape, dtype=np.complex128)
        _hankel1_ratio_2d(ma, x, br)
    
    return br


