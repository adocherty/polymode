from numpy import *
from pylab import *
from Polymode import *

if 'modes' not in locals():
    wl = 1.55
    Nx = (200,21)
    m145 = Material.Fixed(1.45)

    wg = Waveguide.Waveguide(material=m145, symmetry=6)
    wg.add_shape(Waveguide.Circle(Material.Air(), center=(6.75,0), radius=2.5))
    modes = NLSolver.DefaultSolver(wg,Nx)(wl, 0, 1.43, number=2)

#Coordinate object for mode
res = 50
c1 = coordinates.PolarCoord(rrange=(0, 3), border=1, arange=(-pi,pi), N=(res,res))
c2 = coordinates.CartesianCoord(X=10, Y=10, N=res)

m = modes[1]

#The components of e and h
hr,ha,hz = m.magnetic_field(coord=c1)
er,ea,ez = m.electric_field(wg, coord=c1)

print "TM like mode factor", sqrt(abs(hr)**2+abs(ha)**2).sum()/abs(hz).sum()
print "TE like mode factor", sqrt(abs(er)**2+abs(ea)**2).sum()/abs(ez).sum()

#The components of e and h for cartesian fields
h2 = m.magnetic_field(coord=c2)
e2 = m.electric_field(wg, coord=c2)

hx,hy,hz = h2
ex,ey,ez = e2

#The poynting vector
poynting = cross(e2,conj(h2), axis=0)
n2 = wg.index2(wl=wl, coord=m.coord, resample=c2)

#The x-component of the poynting vector
subplot(131)
Plotter.plot_v(c2, poynting[0].real)
wg.plot(fill=0, substrate=0)

subplot(132)
Plotter.plot_v(c2, poynting[1].real)
wg.plot(fill=0, substrate=0)

subplot(133)
Plotter.plot_v(c2, poynting[2].real)
wg.plot(fill=0, substrate=0)
