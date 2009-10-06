from numpy.testing import *

#--------------------------------------------------------------------------
import sys,os,logging
sys.path.append(os.pardir)
#--------------------------------------------------------------------------

from Polymode import *
from Polymode import LayeredSolver

exact_neffs0 = [(1.4631341116878256+0.00078917950990296381j), (1.4629976772638573+0.00080559894572876877j)]
exact_neffs1 = [(1.4743863669498769+0.00016865030027398155j), (1.4486028745001998+0.0021698870899070719j)]

class TesterLayered(TestCase):
    solver_class = LayeredSolver.LayeredSolver

    def _create_problem(self):
        # Materials
        core = Material.SiO2GeO2(0.2)
        cladding = Material.Silica()

        # Create waveguide
        self.wg = Waveguide.Waveguide(material=core, symmetry=1)

        #Create the core
        s1 = Waveguide.Annulus(cladding, r=(2.0,3.0))
        self.wg.add_shape(s1)

        # Create the solver
        self.solver = self.solver_class(self.wg)

    def test_solve(self):
        self._create_problem()
        wl = 1.0
        
        #Solve for 2 modes in mode class 0
        self.modes0 = self.solver(wl, 0, nefflist=[1.46313411166+0.0007891j, 1.46299780529+0.00080562j])

        #Solve for 1 mode in mode class 1
        self.modes1 = self.solver(wl, 1, nefflist=[1.4744+1e-4j, 1.44861+2e-3j])

        neffs0 = [m.neff for m in self.modes0]
        neffs1 = [m.neff for m in self.modes1]

        logging.info("%s" % neffs0)
        logging.info("%s" % neffs1)

        assert_array_almost_equal(neffs0, exact_neffs0)
        assert_array_almost_equal(neffs1, exact_neffs1)

        #Check field profiles

class TesterLayeredCauchy(TesterLayered):
    solver_class = LayeredSolver.LayeredSolverCauchy

if 0:
#class TesterFinitedifference(TesterLayered):
    Nx = (1000,1)
    def _create_problem(self):
        # Materials
        core = Material.SiO2GeO2(0.2)
        cladding = Material.Silica()

        # Create waveguide
        self.wg = Waveguide.Waveguide(material=core, symmetry=1)

        #Create the core
        s1 = Waveguide.Annulus(cladding, r=(2.0,3.0))
        self.wg.add_shape(s1)

        # Create the solver
        self.solver = NLSolver.DefaultSolver(self.wg, self.Nx)


if __name__ == '__main__':
    run_module_suite()

