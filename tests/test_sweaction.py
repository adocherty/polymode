#Unit test for VWE_Action

from __future__ import division

from numpy import *
from numpy.testing import *

#--------------------------------------------------------------------------
#Remove main ABCSolver pathname & add path one up in directory to PATH
import sys,os
#sys.path.append(os.path.split(os.path.abspath(os.curdir))[0])
sys.path.append(os.path.abspath(os.curdir))
#--------------------------------------------------------------------------

from VWE_Action import *

class test_VWE_Action(NumpyTestCase):
	def check_complex(self):
		pass

class test_Hankelratio(NumpyTestCase):
	'''
	Check Hankel ratio function gives correct results
	'''
	def test_hankel(self):
		ms = arange(-200,200,10)
		x = 6.18709499137-0.000430208636459j

		answer = array([-32.30979881 -2.24876725e-03j, -30.69270826 -2.13644030e-03j, -29.07552567 -2.02411976e-03j, -27.45823466 -1.91180680e-03j,-25.84081475 -1.79950283e-03j, -24.22323996 -1.68720969e-03j, -22.60547679 -1.57492971e-03j, -20.98748133 -1.46266598e-03j,-19.36919489 -1.35042262e-03j, -17.75053715 -1.23820530e-03j, -16.13139519 -1.12602199e-03j, -14.51160509 -1.01388429e-03j,-12.89091938 -9.01809806e-04j, -11.26894597 -7.89826599e-04j,  -9.64502393 -6.77982377e-04j,  -8.01794435 -5.66365607e-04j,-6.38522875 -4.55161871e-04j,  -4.74083918 -3.44843658e-04j,  -3.06486731 -2.37145693e-04j,  -1.19981619 +4.75288828e-03j,-0.08032385 +1.00314151e+00j,  -1.19981619 +4.75288828e-03j,  -3.06486731 -2.37145693e-04j,  -4.74083918 -3.44843658e-04j,-6.38522875 -4.55161871e-04j,  -8.01794435 -5.66365607e-04j,  -9.64502393 -6.77982377e-04j, -11.26894597 -7.89826599e-04j,-12.89091938 -9.01809806e-04j, -14.51160509 -1.01388429e-03j, -16.13139519 -1.12602199e-03j, -17.75053715 -1.23820530e-03j,-19.36919489 -1.35042262e-03j, -20.98748133 -1.46266598e-03j, -22.60547679 -1.57492971e-03j, -24.22323996 -1.68720969e-03j,-25.84081475 -1.79950283e-03j, -27.45823466 -1.91180680e-03j, -29.07552567 -2.02411976e-03j, -30.69270826 -2.13644030e-03j])

		self.__check__(ms,x,answer)


	def __check__(self, ms, x, answer):
		h = hankel_ratio(ms, x)
		assert_almost_equal(h, answer)

if __name__ == '__main__':
	NumpyTest().run()
	
