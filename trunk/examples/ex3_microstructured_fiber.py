from numpy import *
from Polymode import *

# Materials
air = Material.Air()
m2 = Material.Fixed(2)
polymer = Material.Polymer()

# Create waveguide
wg = Waveguide.Waveguide(material=polymer, symmetry=6)

# Construct some shapes to put in the waveguide
s1 = Waveguide.Circle(air, center=(3,pi/6), radius=1.25)
wg.add_shapes(s1)

Plotter.figure()
wg.plot()
Plotter.show()
