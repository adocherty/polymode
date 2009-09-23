#Unit test for VWE_Action

from __future__ import division

import numpy as np
from numpy.testing import *

#--------------------------------------------------------------------------
import sys,os
sys.path.append(os.pardir)
#--------------------------------------------------------------------------

from Polymode.mathlink.cauchy_root_solve import *

class test_blockarray(TestCase):
    #1) Delves and Lyness
    polynomial_roots = np.array([1-1j, -0.5-1j, 0.24+0.5j, 0.25+0.5j, 0.245+0.51j, -0.5+0.1j, 0.1+pi*1j, 5+5j, 1+20j])+5.0
    ftest1_roots = np.array([1.41460718-3.04772206j,1.41460718+3.04772206j, -1.84423395,
            0.530894930-1.33179188j, 0.530894930+1.33179188j, 0])

    #Test functions
    def ftest1(self, z):
        f = np.exp(3*z)+2*z*np.cos(z)-1
        return f
    def ftest1prime(self, z):
        fp = 3*np.exp(3*z)+2*np.cos(z)-2*z*np.sin(z)
        return fp

    def fpoly(self, z):
        fx = ones(z.shape, complex)
        for a in self.polynomial_roots:
            fx *= z-a
        return fx
    def fpolyprime(self, z):
        fp = zeros(z.shape, complex)
        fx = self.fpoly(z)

        #As z!=a we can divide
        for a in self.polynomial_roots:
            fp += fx/(z-a)
        return fp

    def Atest(z):
        A = np.eye(2)*z - np.array([[4.5+0.1j,1],[0,4.6-0.1j]])
        return A

    def _compare_best_candidates(self, zs, ans):
        assert_equal(len(zs), len(ans))
        
        sorted_ok = []
        for ii in range(len(ans)):
            aii = ans[ii]
            zarg = np.absolute(zs - aii).argmin()
            sorted_ok.append(zarg)
            
            assert_almost_equal(zs[zarg], aii, decimal=5)

    def _run_algorithms(self, f, fp=None, z0=0, roots=[]):
        if fp is not None:
            print "\nChecking Delves"
            czd = findzero_delves(f, fp, z0=z0, R=4)
            self._compare_best_candidates(czd,roots)

        print "\nChecking  Carpentier"
        czc = findzero_carpentier(f, z0=z0, R=4)
        self._compare_best_candidates(czc,roots)

        print "\nChecking ADR"
        cza = findzero_adr(f, z0=z0, R=4)
        self._compare_best_candidates(cza,roots)

    def test_polynomials(self):
        self._run_algorithms(self.fpoly,self.fpolyprime,5,self.polynomial_roots[:-2])
        
    def test_complexexp(self):
        self._run_algorithms(self.ftest1,self.ftest1prime,0,self.ftest1_roots)

if __name__ == '__main__':
    run_module_suite()
