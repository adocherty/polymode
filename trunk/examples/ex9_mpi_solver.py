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

Nx = 400,15

#Create queue only if we are on the master node
if MPISolver.mpi_is_master():
    #number_of_jobs could be dependent upon
    #the number of processors available
    number_of_jobs = 20
    number_of_modes = 400
    number_of_modes_per_job = number_of_modes//number_of_jobs
    fullneffrange = poly.index(wl), 1.2
    dneff = (fullneffrange[1]-fullneffrange[0])/number_of_jobs

    solvers = []
    for ii in range(number_of_jobs):
        m0 = 0

        #Specify search range
        neffrange = fullneffrange[0] + dneff*ii, fullneffrange[0] + dneff*(ii+1)

        #Create new solver object
        solver = NLSolver.DefaultSolver(wg1, Nx, store=False)
        solver.initialize(wl, m0, neffrange, number=number_of_modes_per_job)

        solvers += [ solver ]

    queue_save_filename = 'mpi_test.queue'

MPISolver.mpi_batch_solve(solvers, queue_save_filename)
