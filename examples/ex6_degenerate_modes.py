from numpy import *
from pylab import *
from Polymode import *

if 'modes' not in locals():
    wl = 1.55
    Nx = (200,21)
    m145 = Material.Fixed(1.45)

    wg = Waveguide.Waveguide(material=m145, symmetry=6)
    wg.add_shape(Waveguide.Circle(Material.Air(), center=(6.75,0), radius=2.5))
    modes = NLSolver.DefaultSolver(wg, Nx)(wl, 1, number=1)

mc1,mc2 = Modes.construct_degenerate_pair(modes[0])

m1,m2 = Modes.construct_lp_degenerate_pair(modes[0])

Plotter.figure()
Plotter.plot_modes_in_grid([mc1,mc2,m1,m2], 'pe')
Plotter.show()
