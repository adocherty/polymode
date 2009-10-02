from pylab import *
from numpy import *
from Polymode import *

air = Material.Air()
m1 = Material.Fixed(1.45)

wg = Waveguide.Waveguide(rmax=9.5, material=m1, symmetry=6)
wg.add_shape( Waveguide.Circle(air, center=(6.75,0), radius=2.5) )

Nx = 200,21                                     #number of radial & azimuthal points
wl = 1.0                                                        #Wavelength

solver = NLSolver.DefaultSolver(wg,Nx)

if 0:
    modes0 = solver(wl, 0, 1.448, number=4)
    modes1 = solver(wl, 1, 1.448, number=4)
    modes2 = solver(wl, 2, 1.448, number=4)
    modes3 = solver(wl, 3, 1.448, number=4)
    modes = modes0+modes1+modes2+modes3

if 1:
    modes3 = solver(wl, 3, 1.448, number=2)
    modesn3 = solver(wl, -3, 1.448, number=2)
    modes = modes3+modesn3

figure(1)
Plotter.plot_modes_in_grid(modes, cmap=cm.hot)
Plotter.plot_modes_in_grid(modes, 'pe')


wgfull = Waveguide.Waveguide(rmax=9.5, material=m1, symmetry=1)
c = Waveguide.Circle(air, center=(6.75,0), radius=2.5)
wgfull.add_shape( Waveguide.create_hexagonal_tiling(c, D=6.75, symmetry=1) )

#Nx = 200,61                                    #number of radial & azimuthal points
#wl = 1.0                                                       #Wavelength

#neffsfull = array([ 1.44441337 +1.25484898e-07j,  1.44432966 +2.32507763e-07j,
#        1.44775232 +6.14743566e-09j,  1.44775232 +6.14743566e-09j,
#        1.44436305 +1.78375392e-07j,  1.44436305 +1.78375392e-07j,
#        1.44033182 +3.61100710e-06j,  1.43997714 +2.07377612e-06j])
#solver = NLSolver.DefaultSolver(wgfull,Nx)
#modesfull = solver(wl, 0, nefflist=neffsfull)

#figure(2)
#Plotter.plot_modes_in_grid(modesfull, cmap=cm.hot)
#Plotter.plot_modes_in_grid(modesfull, 'pe')
