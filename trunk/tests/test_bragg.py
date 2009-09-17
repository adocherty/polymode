from numpy import *
from pylab import *
from mathlink import timer

from PolyMode.LayeredSolver import *

timer.reset()

#Parameters
wl = 1.1
k0 = 2*pi/wl

ng = 1.0
n1 = 1.46
n2 = 1.44

dlayer = 0.5
alayer = (1+sqrt(n1**2-1)/sqrt(n2**2-1))**(-1)
rcore = 5.0
Nlayer = 10
m = 0

#Construct the layers
ns = [ng] + [n1,n2]*Nlayer
gammas = [0]+[0,0]*Nlayer
drs = [rcore] + [alayer*dlayer, (1-alayer)*dlayer]*Nlayer
rs = cumsum(drs[:-1])

#Wavelength Scan
wls = arange(0.8,1.4,1)
Nwl = len(wls)

T = transfer_matrix(m, wl, rs, ns)
neff = find_modes(T, nrange=[0.6,1.0], plot_found='r-')

print timer.report()
