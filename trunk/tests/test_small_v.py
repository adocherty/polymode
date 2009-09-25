import numpy as np
from numpy.testing import *

#--------------------------------------------------------------------------
import sys,os,logging
sys.path.append(os.pardir)
#--------------------------------------------------------------------------

from Polymode import *
from Polymode import LayeredSolver

class tester_layered(TestCase):

    def _create_problem(self, n1, n2, rc):
        # Materials
        core = Material.Fixed(n1)
        cladding = Material.Fixed(n2)

        # Create waveguide
        self.wg = Waveguide.Waveguide(material=cladding, symmetry=1)

        #Create the core
        s1 = Waveguide.Circle(core, radius=rc)
        self.wg.add_shape(s1)

    def test_cauchy(self):
        self._create_problem(1.46,1.44,1.0)

        # Create the solver
        self.solver = LayeredSolver.LayeredSolverCauchy(self.wg)
        self.solver.setup(debug_plot=1, Nscan=10)
        
        #self.solver = LayeredSolver.LayeredSolver(self.wg)
        wls = np.arange(1,3,0.2)
        for wl in wls:
            modes = self.solver(wl, 1, number=1)

    def _test_standard(self):
        self._create_problem(1.46,1.44,1.0)

        # Create the solver
        self.solver = LayeredSolver.LayeredSolver(self.wg)
        self.solver.setup(Nscan=1e3, debug_plot=1)
        
        #self.solver = LayeredSolver.LayeredSolver(self.wg)
        wls = np.arange(1,3,0.2)
        for wl in wls:
            modes = self.solver(wl, 1, number=1)




if __name__ == '__main__':
    run_module_suite()

