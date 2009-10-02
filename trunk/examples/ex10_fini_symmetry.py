from __future__ import division
from pylab import *
from numpy import *

from Polymode import *
from Polymode import KrylovSolver

## Waveguide Parameters:
Nx=300,20
wl = 1.44

## Materials
air = Material.Air()
polymer = Material.Polymer()

## Create waveguide
wg = Waveguide.Waveguide(rmin=0, rmax=12, material=polymer, symmetry=6)

#Construct the waveguide
core = Waveguide.Circle(air, center=(0,0), radius=3.0)
c = Waveguide.Circle(air, radius=1.0)
shapes = Waveguide.create_hexagonal_tiling(c, offset=2, layers=3, D=2.5, symmetry=6)

wg.add_shapes(core, shapes)

solver = NLSolver.DefaultSolver(wg, Nx)

if 0:
    mds=[]
    wls = arange(1.4,1.48,0.001)
    for wl in wls:
        mds += solver(wl, 0, 0.98, number=10)
        mds += solver(wl, 1, 0.98, number=10)
        mds += solver(wl, 2, 0.98, number=10)

elif 1:
    #Create WL range solver
    wlsolver = Solver.WavelengthTrackingSolver(solver)
    wlsolver.ga_target = 1e-4

    nefflist = [0.9999, 0.983428, 0.986677, 0.97, 0.965, 0.96]
    wlrange=[1.4,1.55]
    modes0 = wlsolver(wlrange, 0, nefflist=nefflist)
    #modes1 = wlsolver(wlrange, 1, nefflist=[0.99, 0.995])

elif 1:
    wg = Waveguide.Waveguide(rmin=2.0, rmax=12, material=polymer, symmetry=6)
    wg.add_shapes(core, shapes)

    #Create WL range solver
    ksolver = KrylovSolver.NLBiArnoldiSolver(wg, Nx)
    wlsolver = Solver.WavelengthTrackingSolver(ksolver)
    wlsolver.ga_target = 1e-4
    nefflo = [0.9999, 0.983428, 0.986677, 0.97, 0.965, 0.96]

    modes=[]
    wlrange=[1.4,1.55]
    for neffx in nefflo:
        modes += wlsolver(wlrange, 0, nefflist=[neffx])

else:
    from PolyMode.ExternalSolvers import *

    jcmsolver = JCMWaveSolver(wg, numthreads=4, project_path='jcm', \
            main_path='/home/andrewd/Software/JCMsuite')
    jcmsolver.setup(xy_domain=1, geo_accuracy=2, refinement_levels=3)
    jcmsolver.export_triangulation(view=1)
