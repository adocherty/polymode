from numpy import *
from Polymode import *

## Materials
air = Material.Air()
m2 = Material.Fixed(2)
polymer = Material.Polymer()

## Create waveguide
wg = Waveguide.Waveguide(material=polymer, symmetry=3)

##Construct circle in front of ellipse using zorder
s0 = Waveguide.Ellipse(m2, center=(0,0), axes=(3.0,1.5))
s1 = Waveguide.Circle(air, center=(3,pi/6), radius=1.5, zorder=1)
wg.add_shapes(s0, s1)

Plotter.figure()

Plotter.subplot(121)
wg.plot()

#Change orders around
s0.zorder=1
s1.zorder=0

Plotter.subplot(122)
wg.plot()

Plotter.show()
