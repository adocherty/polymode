import logging, datetime
#logging.basicConfig(level=logging.DEBUG, format='%(created)s|%(process)d|%(levelname)s|%(message)s', filename="mpi-%s.log" % datetime.date.today() )

from numpy import *
from Polymode import MPISolver, Material, Waveguide, NLSolver

#Materials
air = Material.Air()
poly = Material.Polymer()

#Create a waveguide
symmetry = 10
wl = 1.0
d = 1.0
w = 5.0
rcore = 40
rmax = 45
wedge1 = Waveguide.ChannelAnnulus(air, r=(rcore,rcore+w), phi=(-pi/symmetry,pi/symmetry), d=d)

wg1 = Waveguide.Waveguide(material=poly, symmetry=symmetry)
wg1.add_shape(wedge1)

Nx = 100,11
#Create solvers

#number_of_jobs could be dependent upon
#the number of processors available
if 1:
    number_of_jobs = 5
    totalnumber = 20
    fullneffrange = poly.index(wl), 1.2
    dneff = (fullneffrange[1]-fullneffrange[0])/number_of_jobs

    solvers = []
    for ii in range(number_of_jobs):
        m0=0

        #Specify search range
        neffrange = fullneffrange[0] + dneff*ii, fullneffrange[0] + dneff*(ii+1)

        #Create new solver object
        solver = NLSolver.DefaultSolver(wg1, Nx)
        solver.initialize(wl, m0, neffrange, number=totalnumber//number_of_jobs)

        solvers += [ solver ]

    queue_save_filename = 'mpi_test.queue'
    MPISolver.mpi_batch_solve(solvers, queue_save_filename)

