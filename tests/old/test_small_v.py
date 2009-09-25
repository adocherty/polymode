from numpy import *
from pylab import *

from PolyMode.LayeredSolver import *

#Parameters
m0 = 1
wl = 1.0
k0=2*pi/wl

n1 = 1.45
n2 = 1.44
rs = arange(0.4,1.5,0.1)
for rc in rs:
	T = transfer_matrix(m0, wl, [rc], [n1,n2])
	v = 2*pi*rc*sqrt(n1**2-n2**2)/wl
	print "Finding modes for V=%.2g" % v

	neff = find_modes(T, res=1e3, plot_found='k:')
	
pl.draw()
pl.show()
