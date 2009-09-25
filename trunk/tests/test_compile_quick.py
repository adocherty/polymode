
from numpy.testing import *

#--------------------------------------------------------------------------
import sys,os
sys.path.append(os.pardir)
#--------------------------------------------------------------------------

from Polymode.mathlink.cauchy_root_solve import *

class test_compile(TestCase):
    def test_bessel_compile(self):
        from Polymode.mathlink import bessel_ratios

    def test_ublocklu_compile(self):
        from Polymode.mathlink import ublocklu

if __name__ == '__main__':
    run_module_suite()
    
