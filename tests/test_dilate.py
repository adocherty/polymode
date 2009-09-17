#test_waveguide.py

from pylab import *

from ABCSolver import *
from ABCSolver.Waveguide import *

air = Material.Air()
m1 = Material.Silica()
m2 = Material.Gold()

Nx = 100,41
wg = Waveguide(material=air, symmetry=2)
coord = wg.get_coord(Nx)
ac = AnnularCombine(coord,2)

c1 = Circle(m1, center=(2.5,pi/3), radius=0.5)
c2 = Circle(m2, center=(2.5,pi/3), radius=0.5)
c2.zorder=-1
c2.dilate(0.1)
wg.add_shape( c1, c2 )

#Note: Nodes must be arranged in clockwise order for dilate to work!
polyg = Polygon(m1, nodes=[(1,0),(2,0),(1.5,pi/6)])
polyg2 = Polygon(m2, nodes=[(1,0),(2,0),(1.5,pi/6)])
polyg2.zorder = -1
polyg2.dilate(0.1)
wg.add_shape( polyg, polyg2 )

rect = Rectangle(m1, center=(1,pi/3), axes=(0.6,0.6))
rect2 = Rectangle(m2, center=(1,pi/3), axes=(0.6,0.6))
rect2.zorder = -1
rect2.dilate(0.1)
wg.add_shape( rect, rect2 )

ell = Ellipse(m1, center=(2,pi/3), axes=(0.5,0.25), rot=0, xy=1)
ell2 = Ellipse(m2, center=(2,pi/3), axes=(0.5,0.25), rot=0, xy=1)
ell2.zorder = -1
ell2.dilate(0.1)
wg.add_shape( ell, ell2 )

ann = Annulus(m1, r=(1.5,2.5), phi=(-pi/6, -pi/3))
ann2 = Annulus(m2, r=(1.5,2.5), phi=(-pi/6, -pi/3))
ann2.dilate(0.2)
ann2.zorder=-1
wg.add_shape( ann, ann2 )

clf()
wg.plot()
show()
