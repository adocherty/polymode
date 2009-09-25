from __future__ import division

from pylab import *
from numpy import *
from numpy import testing

from ABCSolver import *

class TestDiff(testing.TestCase):

	def test_image_load(self):
		air = Material.Air()
		m1 = Material.Fixed(1.5)

		wl=1.0
		Nx = 100,21
		icenter = array([ 342.34780297,  345.38735766])

		#Test waveguide
		wg_test = Waveguide.Waveguide(material=air, symmetry=6)
		coord = wg_test.get_coord(Nx)

		image = Image.ImportBitmap("MOF_NRL.png", size=(60,55), negative=False)
		image.contrast(scale=5, transition=0.6, clip=1)

		#Check center calculation
		center = image.get_center(coord)
		testing.assert_array_almost_equal(icenter, center)
		
		#Create waveguide shape
		wg = Waveguide.Waveguide(rmax=5, material=air, symmetry=6)
		wg.add_shape(Waveguide.Image(m1, image))
		inx = wg.index(wl, Nshape=(10,11))
		inx_true = array([[ 1.47974688,  1.47974688,  1.47974688,  1.47974688,  1.47974688,  1.47974688,  1.47974688,  1.47974688,  1.47974688,  1.47974688,  1.47974688],
       [ 1.46013204,  1.46735661,  1.46300012,  1.46865651,  1.48668883,  1.48869629,  1.49174347,  1.49174347,  1.49174347,  1.47464849,  1.48344089],
       [ 1.46739475,  1.42027624,  1.49578973,  1.5       ,  1.5       ,  1.46060734,  1.49174347,  1.49174347,  1.49174347,  1.49174347,  1.5       ],
       [ 1.36428992,  1.28668098,  1.04970456,  1.05652995,  1.0091836 ,  1.01050747,  1.33401798,  1.14354776,  1.43053026,  1.5       ,  1.48096587],
       [ 1.08534806,  1.06589942,  1.02931823,  1.01055044,  1.00277665,  1.00091266,  1.00000125,  1.00000208,  1.        ,  1.00003375,  1.12990853],
       [ 1.11999304,  1.06573207,  1.01945692,  1.00269217,  1.00111502,  1.00011819,  1.00000625,  1.0000037 ,  1.00000216,  1.00000612,  1.08592044],
       [ 1.43239834,  1.07594597,  1.02521715,  1.00317637,  1.00054492,  1.00008216,  1.00006131,  1.00001221,  1.0000038 ,  1.00000698,  1.00022721],
       [ 1.5       ,  1.5       ,  1.03613905,  1.00196067,  1.00029564,  1.00003739,  1.00097994,  1.00007696,  1.00002801,  1.00000934,  1.5       ],
       [ 1.0000259 ,  1.00002035,  1.5       ,  1.03873628,  1.00589337,  1.00020212,  1.00366748,  1.00067495,  1.00009212,  1.19118449,  1.00225351],
       [ 1.00003116,  1.00000327,  1.00000386,  1.5       ,  1.11940671,  1.02662355,  1.5       ,  1.5       ,  1.5       ,  1.03051892,  1.00093383],
       [ 1.00007016,  1.000006  ,  1.00000402,  1.29691924,  1.03673579,  1.00138188,  1.00000041,  1.00000003,  1.5       ,  1.04220741,  1.00231671],
       [ 1.00009757,  1.00000888,  1.00150987,  1.38073718,  1.01082597,  1.00071852,  1.00002455,  1.00000005,  1.00044301,  1.07454852,  1.00370667]])

		testing.assert_array_almost_equal(inx, inx_true)

		#Find mode
		wg = Waveguide.Waveguide(rmax=10, material=air, symmetry=6)
		wg.add_shape(Waveguide.Image(m1, image))
		bs = Solver.DefaultSolver(Nx, wg)
		md, = bs(1, wl, 1.47, number=1)
		neff = 1.4710337210660027+1.8972754998663313e-05j
		testing.assert_almost_equal(md.neff, neff)


if __name__ == "__main__":
		testing.run_module_suite()

