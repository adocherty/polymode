#Unit test for hankel
from __future__ import division
import sys,os

from PolyMode import mathlink
reload(mathlink)

import pylab as pl
from numpy import *
from numpy.testing import *

from scipy.special import *
from PolyMode.mathlink import timer, besselj_ratio, hankel1_ratio

x = 83+1j
m = arange(-300,300)

tick = timer.timer()

#Test bessel ratio functions
tick.start("scipy")
br1 = x*jvp(m,x)/jv(m,x)
tick.stop()

tick.start("br")
br2 = besselj_ratio(m,x)
tick.stop()

print "Max error in bessel ratio", abs(br1-br2).max()

tick.start("scipy")
hr1 = m - x*hankel1(m+1,x)/hankel1(m,x)
tick.stop()

tick.start("hr")
hr2 = hankel1_ratio(m,x)
tick.stop()

print "Max error in hankel ratio", abs(hr1-hr2).max()

#pl.plot(m, abs(hr1-hr2), 'r-')

tick.report()
