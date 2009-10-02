from numpy import *
from Polymode import Material, Plotter

#To search the datafiles for a specific material  use the find_materials function
Material.find_materials('ag')

#Use supplied data from internal material databases
ag_fs = Material.FreeSnell('ag')

#Compare refractive index of some plastics
pmma = Material.PMMA()
pc = Material.Polycarbonate()
ps = Material.Polystyrene()
zeo = Material.Zeonex()

#Plot the dispersion of the different materials
#This will give some errors as we are extrapolating
#outside the valid region of wavelengths
for mat in [pmma, pc, ps, zeo]:
    print mat
    mat.plot([0.4,1.2], showdata=True)
Plotter.title("The refractive index of different optical polymers")
Plotter.show()

#Example of a new material using Sellmeier coefficients
B = [0.6961663, 0.4079426, 0.8974794]
L = [0.0684043, 0.1162414, 9.896161]
nm_sellmeier = Material.Sellmeier(B,L)

#Example of a new material using a Laurent series
A = [-1.25900000e-02, 2.38490000, 1.07900000e-02, 1.65180000e-04, \
                -1.94741000e-06, 9.36476000e-08 ]
nm_laurant = Material.Laurent(A)

#Example of a new custom material using a formula
def strange_index(self, wavelength):
    return 1+0.1*wavelength**2

nm_formula = Material.Material()
nm_formula.index_function = strange_index
nm_formula.name = "Strange Function Material"
nm_formula.color = 'blue'

#Example of a new material using interpolation of a Sopra formatted file
#nm_sopra = Material.SopraFile( open('test.nk', 'r') )
