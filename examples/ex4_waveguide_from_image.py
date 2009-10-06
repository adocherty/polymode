from numpy import *
from Polymode import *

import os

# Materials
air = Material.Air()
m2 = Material.Fixed(2)
polymer = Material.Polymer()

image_filename1 = os.path.join("images","MOF_NRL.png")
image_filename2 = os.path.join("images","hlfcollapsedKagome.png")

# Create waveguide and fix the radial limit to be 15
wg = Waveguide.Waveguide(rmax=5, material=polymer, symmetry=6)

## Load a png image to put into the waveguide
image = Image.ImportBitmap(image_filename1, size=(20,20), offset = pi/6, negative=True)
image.contrast(scale=15, transition=0.45)
wg.add_shape(Waveguide.Image(air, image))
print wg

# Another example of bitmap import
wg2 = Waveguide.Waveguide(material=air, symmetry=6)
kagim = Image.ImportBitmap(image_filename2, size=(100,100), negative=True)
wg2.add_shape(Waveguide.Image(polymer, kagim, average_over_sectors=1))
print wg2

Plotter.figure()

# Plot two waveguides
Plotter.subplot(121)
wg.plot_bitmap((200,60), sectors=1, cmap=Plotter.cm.gray)

Plotter.subplot(122)
wg2.plot_bitmap((400,100), sectors=1, cmap=Plotter.cm.gray)

Plotter.show()
