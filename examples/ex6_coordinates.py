from Polymode import *

from numpy import *
from pylab import *

res = 20

#Coordinate object for mode
c1 = coordinates.PolarCoord(rrange=(0.0,2), border=1, arange=(-pi,pi), N=(res,2*res))
c2 = coordinates.CartesianCoord(X=2, Y=2, N=res)

#Get r, phi
r1,phi1 = c1.polar2d()
r2,phi2 = c2.polar2d()

#Get cartesian points
x1,y1 = c1.cartesian2d()
x2,y2 = c2.cartesian2d()

#Plot the point distribution
subplot(221)
scatter(x1.ravel(), y1.ravel(), 2)
axis('tight')

subplot(222)
scatter(x2.ravel(), y2.ravel(), 2)
axis('tight')

#Example function
f = lambda x,y: y*exp(-2*(x**2 + y**2))

g1=c1.div_t(c1.grad_t(f(x1,y1)))
g2=c2.div_t(c2.grad_t(f(x2,y2)))

subplot(223)
Plotter.plot_v(c1,real(g1))
axis('tight')

subplot(224)
Plotter.plot_v(c2,real(g2))
axis('tight')

print "Function integral on polar grid:", c1.int_dA(abs(f(x1,y1)))
print "Function integral on cartesian grid:", c2.int_dA(abs(f(x2,y2)))

draw()
