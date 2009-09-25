# _*_ coding=utf-8 _*_
import logging
from numpy import *

from scipy import fftpack as fftp

_output_unicode = False
def utf8out(s):
    """
    Encode string to utf8 if unicode output is enabled, otherwise replace with ascii
    representations of technical symbols.
    """
    if _output_unicode:
        s = s.encode('utf8')
    else:
        unicode_greek = {'α':'alpha', 'β':'beta','γ':'gamma','δ':'delta','ε':'epsilon',
            'θ':'theta', 'λ':'wl', 'μ':'u', 'π':'pi', 'φ':'phi', 'ω':'ang. freq.', 'ϕ':'phi',
            '₀':'0', '₂':'2', '₁':'1', '₃':'3', '₄':'4', '₅':'5', '₆':'6', '₇':'7', '₈':'8', '₉':'9'}
        
        for x in unicode_greek.keys():
            s = s.replace(unicode(x,'utf8'), unicode_greek[x])

        s = s.encode('ascii', 'ignore')
    return s

#Import bessel functions from scipy or amos
#from .amos import hankel1, hankel1p, hankel2, hankel2p

from scipy.special import hankel1, hankel2, jv, jvp, kv, kvp
def hankel1p(m,z):
    return hankel1(m-1,z)-hankel1(m,z)*m/z
def hankel2p(m,z):
    return hankel2(m-1,z)-hankel2(m,z)*m/z

#Accurate Bessel function ratio functions
try:
    from .bessel_ratios import besselj_ratio, hankel1_ratio

except ImportError:
    logging.warning("Using scipy hankel ratio - this may be inaccurate!")
    def hankel1_ratio(m, x):
        "Calculates the ratio of Hankel functions: hankel_ratio = hankel1'(m,x)/hankel1(m,x)"
        br = x*hankel1(m-1,x)/hankel1(m,x) - m
        if any(isnan(br)) or any(isinf(br)): logging.error("Nan in Scipy Bessel: compile bessel module")
        return br

    #The DtN function for the inner boundary condition
    def besselj_ratio(m, x):
        "Calculates the ratio of Hankel functions: hankel_ratio = hankel1'(m,x)/hankel1(m,x)"
        br = x*jv(m-1,x)/jv(m,x) - m
        if any(isnan(br)) or any(isinf(br)): logging.error("Nan in Scipy Bessel: compile bessel module")
        return br

#Import from scipy otherwise fallback to internal
#All choices should be put here, if scipy is found and defaults if not
try:
    from scipy import constants
except ImportError:
    from . import physcon as constants

__all__ = filter(lambda s:not s.startswith('_'),dir())


