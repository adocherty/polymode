import os
from numpy import *
from Polymode import *

# Materials
air = Material.Air()
polymer = Material.Fixed(1.49)

image_file_name = os.path.join("images","MOF_NRL.png")

# Create waveguide and specify the problem domain to have radius 15
# With an Image shape we need to specify the interior material
#wg = Waveguide.Waveguide(rmax=15, material=polymer, symmetry=1)

# Load a png file and specify the size of the original image
# The center of the image (pixel offset) must also be specified for symmetry=1
#imagefilename = os.path.join("images","celery.png")
#image = Image.ImportBitmap(imagefilename, size=(60,120), center=(230,140), negative=True)

wg = Waveguide.Waveguide(rmax=5, material=polymer, symmetry=6)
image = Image.ImportBitmap(image_file_name, size=(40,40), offset = pi/6, negative=True)
image.polar_regrid = Image.nn_polar_regrid
image.contrast(scale=15, transition=0.45)

# Add the shape to the waveguide
wg.add_shape(Waveguide.Image(air, image))

# Waveguide Parameters:
Nx=100,60                       #number of radial & azimuthal points
wl=1.55                         #Wavelength

# Create the solver
solver = NLSolver.DefaultSolver(wg, Nx)
modes = solver(wl, 1, 1.4899, number=2)
modes += solver(wl, 0, 1.4899, number=2)

Plotter.figure()
Plotter.subplot(121)
wg.plot_bitmap((100,30), style='lines', alpha=0.5, cmap=Plotter.cm.gray)
modes[0].plot(cmap=Plotter.cm.hot)

Plotter.subplot(122)
wg.plot_bitmap((100,30), style='lines', alpha=0.5, cmap=Plotter.cm.gray)
modes[1].plot(cmap=Plotter.cm.hot)

Plotter.show()
