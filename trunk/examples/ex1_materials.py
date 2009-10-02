from numpy import *
from Polymode import Material, Plotter

#Basic materials
air = Material.Air()
m145 = Material.Fixed(1.45)

#Substrate materials
silica = Material.Silica()
germania = Material.Germania()
pmma = Material.Polymer()

# Accurate ref index for silica doped with 5% Germanium
sige = Material.SiO2GeO2(0.05)

wl = 1.45       #Wavelength in um

for mat in [pmma, silica, sige]:
    ri = mat.index(wl)

    #Get some info
    mat.info()

    #Plot material refractive index
    mat.plot()

    #calculate the refractive index
    print "Refractive index at %gum is %.6g+%.4gi\n" % (wl, real(ri), imag(ri))
