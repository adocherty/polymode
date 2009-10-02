from pylab import *
from numpy import *
from Polymode import *

#Only solve for modes if we don't have them
#Try running this file in ipython twice, with the -i
#option - the second time the modes already calculated
#will be used
if 'modes' not in locals():
    wl = 1.55
    Nx = (400,41)
    m145 = Material.Fixed(1.45)

    wg = Waveguide.Waveguide(material=m145, symmetry=6)
    wg.add_shape(Waveguide.Circle(Material.Air(), center=(6.75,0), radius=2.5))

    solver = NLSolver.DefaultSolver(wg, Nx)
    modes = solver(wl, 1, number=1)

#Coordinate object for mode
res = 20
c1 = coordinates.PolarCoord(rrange=(0,10), border=1, arange=(-pi,pi), N=(res,2*res))
c2 = coordinates.CartesianCoord(X=5, Y=5, N=res)

#Select the first mode
m = modes[0]

#Plot the resampled modal fields
figure()
subplot(121)
m.plot(coord=c1)
subplot(122)
m.plot(coord=c2)

#Mode quantities can be calculated with a user-defined coordinate
print "The group index using polar coordinates:", m.group_index(wg, coord=c1).real
print "The group index using cartesian coordinates:", m.group_index(wg, coord=c2).real

#The transverse components of the electric and magnetic fields
ht = m.magnetic_transverse_field(coord=c1)
et = m.electric_transverse_field(coord=c1)
exh = c1.cross(et, conj(ht))

#The electric and magnetic fields
h = m.magnetic_field(coord=c1)
e = m.electric_field(wg, coord=c1)
h2 = sum(abs(h)**2,axis=0)
e2 = sum(abs(e)**2,axis=0)

#The refractive index squared
n2 = m145.index(wl)**2

#The integral approximation to the propagation constat
beta = 2*(2*pi/wl)*c1.int_dA(n2*exh)/c1.int_dA(h2 + conj(n2)*e2)

print "Propagation constant from integral", beta
print "Propagation constant of mode", m.beta
