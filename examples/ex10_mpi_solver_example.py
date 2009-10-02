"""
MPI Solver example

"""
import logging, datetime
logging.basicConfig(level=logging.INFO, format='%(created)s|%(process)d|%(levelname)s|%(message)s', filename="mpi-%s.log" % datetime.date.today() )

from numpy import *
from Polymode import MPISolver, Material, Waveguide, Solver

#Materials
air = Material.Air()
poly = Material.Polymer()

#Create a waveguide
symmetry = 20
d = 1.0
w = 1
rcore = 40
rmax = 45
wedge1 = Waveguide.ChannelAnnulus(air, r=(rcore,rcore+w), phi=(-pi/symmetry,pi/symmetry), d=d)

wg1 = Waveguide.Waveguide(rmax, poly, symmetry=symmetry )
wg1.add_shape(wedge1)

#Create solvers
Nx = 100,11
k0 = 2*pi/1.0

#number_of_jobs could be dependent upon
#the number of processors available
number_of_jobs = 20
totalnumber = 500
fullneffrange = poly.index(k0), 1.0
dneff = (fullneffrange[1]-fullneffrange[0])/number_of_jobs

solvers = []
for ii in range(number_of_jobs):
    m0=0

    #Specify search range
    neffrange = fullneffrange[0] + dneff*ii, fullneffrange[0] + dneff*(ii+1)

    #Create new solver object
    solver = Solver.DefaultSolver(Nx, wg1)
    solver.initialize(m0, k0, neffrange, number=totalnumber//number_of_jobs)

    solvers += [ solver ]

queue_save_filename = 'data/mpi_test.queue'
MPISolver.mpi_batch_solve(solvers, queue_save_filename)
