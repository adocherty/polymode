from numpy import *
from Polymode import *

# Materials
core = Material.SiO2GeO2(0.2)
cladding = Material.Silica()

# Create waveguide
wg = Waveguide.Waveguide(material=cladding, symmetry=1)

#Create the core
s1 = Waveguide.Circle(core, center=(0,0), radius=2.0)
wg.add_shape(s1)

#Plot the waveguide
Plotter.figure()
wg.plot()
Plotter.show()
