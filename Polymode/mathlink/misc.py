# _*_ coding=utf-8 _*_
'''
Miscellaneous mathematical functions
'''
import numpy as np

#Miscellaneous functions
def terminal_color(msg, color=None):
    if color is None:
        color = 'black'
    colors = {'black':30,'red':31, 'green':32, 'yellow':33, 'blue':34, 'magenta':35,
            'cyan':36, 'white':37, 'reset':39}
    color_str = lambda x,s: "\033[" + str(s) + ";" + str(colors[x.lower()]) + "m"
    return color_str(color,1) + msg + color_str('reset',0)

def format_complex(x, pres=(6,2)):
    if np.isnan(x): return "nan"
    fs = "%%.%dg %%s %%.%dei" % pres
    signs = {1:"+", 0:"", -1:"-"}[np.sign(np.imag(x))]
    return fs % (np.real(x),signs,abs(np.imag(x)))

def format_complex_latex(x, pres=(6,2)):
    if np.isnan(x): return "nan"
    
    fs = r"%%.%dg %%s %%.%df\times 10^{%%d}i" % pres
    signs = {1:"+", 0:"", -1:"-"}[np.sign(np.imag(x))]
    imagexp = np.floor(np.log10(abs(np.imag(x))))
    imagcoef = abs(np.imag(x))/10**imagexp
    
    return fs % (np.real(x),signs,imagcoef,imagexp)

def numands(n):
    sdict = {False:"s", True:""}
    return (n,sdict[n==1])

## Miscellaneous mathematical functions

def sech(x):
    "Returns the hyperbolic secant of x"
    return 1./np.cosh(x)

def coth(x):
    return 1./np.tanh(x)

def absmax(*args):
    if len(args)==1:
        if np.iterable(args[0]):
            args = np.ravel(args[0])
    return args[np.absolute(args).argmax()]

def absmin(*args):
    if len(args)==1:
        if np.iterable(args[0]):
            args = np.ravel(args[0])
    return args[np.absolute(args).argmin()]

def asign(f):
    '''Return the signed magnitude of f. If f is purely real or imaginary then
    this is equivalent to real(f) or imag(f) respectively. Otherwise f is the
    absolute value multiplied by the sign of the largest component.'''
    if np.absolute(np.real(f)).sum() > np.absolute(np.imag(f)).sum():
        af = np.absolute(f)*np.sign(np.real(f))
    else:
        af = np.absolute(f)*np.sign(np.imag(f))
    return af

def machine_precision(dtype=float):
    "Caclulated the machine precision of a number with type dtype"
    d=dtype(1.0)
    while (1.0+d)!=1.0:
        d/=2
    return d

def bracket(x):
    return minimum(1.0,maximum(0.0,x))

def rotationpmatrix(rot):
    rotmat = array([[cos(rot),-sin(rot)],[sin(rot),cos(rot)]])
    return rotmat

def shift_matrix(x, row=False, direction=1, initial=0, dtype=float):
    N = len(x)
    A = np.empty((N,N), dtype=dtype)

    for j in np.arange(0,N):
        if row:
                A[j] = np.append(x[-direction*j+initial:], x[:-direction*j+initial])
        else:   
                A[:,j] = np.append(x[-direction*j+initial:], x[:-direction*j+initial])
    return A
