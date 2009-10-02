from pylab import *
from numpy import *
from Polymode import *

## Materials
air = Material.Air()
m2 = Material.Fixed(2.0)
polymer = Material.Polymer()

## Create waveguide
wg = Waveguide.Waveguide(material=polymer, symmetry=6)

#Construct some shapes to put in the waveguide
s1 = Waveguide.Ellipse(air, center=(3,pi/6), axes=(1.5,1.25))
s2 = Waveguide.Ellipse(air, center=(6,0), axes=(1.5,1.25))
s3 = Waveguide.Annulus(m2, r=(5,7), phi=(pi/12,pi/3-pi/12))

wg.add_shape(s1, s2, s3)

#Create the solvers in a list
wl = 1.45
Naz = 41
solvers = []
Nrs = [50,100,200,400,800]
for Nr in Nrs:
    solver = NLSolver.DefaultSolver(wg, (Nr,Naz))
    solver.initialize(wl, 1, 1.45, number=1)
    solvers.append(solver)

modes = Solver.batch_solve(solvers, filename="convergencerun.solve")

#Extract neff and Nr
neff = [m.neff for m in modes]

figure()

#Plot Nr vs convergence
plot(Nrs[:-1], abs(neff [:-1] - neff[-1]), 'ro')
loglog()
loglog()
xlabel(r'$N_r$')
ylabel(r'Error in $n_{eff}$')

show()
