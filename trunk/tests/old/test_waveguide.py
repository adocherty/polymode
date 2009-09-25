# Test waveguide functions
#
# Todo:
# * Check more shapes
# * Check drlogn and dalogn
#
from numpy import *
from numpy import testing
from ABCSolver import *
from ABCSolver.Waveguide import *

class TestDiff(testing.TestCase):

	def _test_index(self):
		m1 = Material.Fixed(1.0)
		m2 = Material.Fixed(2.0)
		m3 = Material.Fixed(3.0)

		circle = Circle(m2, center=(1.5,0), radius=0.75)
		polyg = Polygon(m2, nodes=[(1,0),(2,0),(1.5,pi/6)])
		rect = Rectangle(m2, center=(2,pi/3), axes=(1.2,1.2))
		ellipse = Ellipse(m2, center=(2,pi/3), axes=(0.5,0.25), rot=-pi/6, zorder=1)
		annulus = Annulus(m2, r=(2.2,3.), phi=(-pi/6, pi/3), zorder=2)

		wg = Waveguide(material=m1, symmetry=3)
		wg.shapes = [circle]
		n2 = wg.index2(1.0, Nshape=(10,11))
		n2_true =array([[ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 1.542434,  2.029663,  1.542434,  0.992134,  1.006509,  0.99689 ,  1.000922,  1.000922,  0.99689 ,  1.006509,  0.992134],
			   [ 4.200297,  3.837064,  4.200297,  2.952429,  0.975401,  1.022645,  0.993162,  0.993162,  1.022645,  0.975401,  2.952429],
			   [ 3.865529,  4.12523 ,  3.865529,  4.156258,  1.626661,  0.888257,  1.026221,  1.026221,  0.888257,  1.626661,  4.156258],
			   [ 3.835908,  4.150785,  3.835908,  4.215322,  1.700764,  0.874764,  1.029477,  1.029477,  0.874764,  1.700764,  4.215322],
			   [ 4.098432,  3.926072,  4.098432,  3.646715,  1.117072,  0.986318,  1.002348,  1.002348,  0.986318,  1.117072,  3.646715],
			   [ 4.312576,  3.7458  ,  4.312576,  2.42749 ,  0.853211,  1.054849,  0.984959,  0.984959,  1.054849,  0.853211,  2.42749 ],
			   [ 2.828561,  4.15977 ,  2.828561,  1.032495,  1.00777 ,  0.99512 ,  1.001536,  1.001536,  0.99512 ,  1.00777 ,  1.032495],
			   [ 1.000459,  1.022136,  1.000459,  0.999894,  1.000042,  0.999982,  1.000005,  1.000005,  0.999982,  1.000042,  0.999894],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ]] , dtype=complex_)

		testing.assert_array_almost_equal(n2, n2_true)

		wg = Waveguide(material=m1, symmetry=3)
		wg.shapes = [polyg]
		n2 = wg.index2(1.0, Nshape=(10,11))
		n2_true = array([[ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 0.777574,  2.036716,  1.917269,  0.801023,  1.129213,  0.89834 ,  1.089803,  0.913266,  1.091022,  0.895276,  1.136423],
			   [ 0.793543,  2.41426 ,  4.385881,  2.057672,  0.942566,  1.002951,  1.016309,  0.971715,  1.040774,  0.94031 ,  1.097568],
			   [ 0.704664,  2.522128,  4.2575  ,  3.688565,  1.609351,  0.798639,  1.142127,  0.875623,  1.123702,  0.861986,  1.177681],
			   [ 0.662871,  2.62002 ,  3.380507,  0.919928,  1.144159,  0.873409,  1.117458,  0.882892,  1.126093,  0.851304,  1.199183],
			   [ 0.848943,  1.891186,  1.315721,  0.878166,  1.080404,  0.936038,  1.056968,  0.94458 ,  1.058598,  0.931957,  1.08989 ],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ],
			   [ 1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ,  1.      ]], dtype=complex_)

		testing.assert_array_almost_equal(n2, n2_true)


	def index_resample(self, wg, Nx):
		c1 = wg.get_coord(Nshape=Nx)
		c2 = coordinates.PolarCoord(rrange=(0,wg.rmax+2*c1.dr),
			arange=(-pi/wg.symmetry,pi/wg.symmetry), N=(Nx[0]+2,Nx[1]), border=0)

		#Get raw index and compare to resampled index
		inx1 = wg.index(1.0, coord=c1)
		inx2 = wg.index(1.0, coord=c1, resample=c2)

		#Note that the raw index returned by wg.index is in the DFT shifted format
		testing.assert_array_almost_equal(fft.fftshift(inx1, axes=[1]), inx2)
		#testing.assert_array_almost_equal(inx1, inx2)

	def test_resample(self):
		m1 = Material.Air()
		m2 = Material.Fixed(1.5)

		for Nx in [(50,21), (50,20), (49,16)]:
			for symm in [10,11]:
				wg = Waveguide(rmin=0, material=m1, symmetry=symm)
				wg.add_shape(Polygon(m2, nodes=[(2,-0.2),(3,-0.4),(3.2,0.5)], xy=1))
				self.index_resample(wg, Nx)
				
				wg = Waveguide(rmin=0, material=m1, symmetry=symm)
				wg.add_shape(Circle(m2, center=(3.0, 0.15), radius=0.25))
				self.index_resample(wg, Nx)

#Manual testing
if 1:
	import pylab as pl
	
	m1 = Material.Air()
	m2 = Material.Fixed(1.5)
	Nx=(50,11)
	wg = Waveguide(rmin=0, material=m1, symmetry=10)
	wg.add_shape(Polygon(m2, nodes=[(2,-0.2),(3,-0.4),(3.2,0.5)], xy=1))
				
	c1 = wg.get_coord(Nshape=Nx)
	rm1,phim1 = c1.polar2d()
	inx1 = fft.fftshift(wg.index(1.0, coord=c1), axes=[1])

	c2 = coordinates.PolarCoord(rrange=(0,wg.rmax+2*c1.dr), arange=(-pi/10,pi/10), N=(52,21), border=0)
	rm2,phim2 = c2.polar2d()
	inx2 = wg.index(1.0, coord=c1, resample=c2)

	pl.subplot(211)
	pl.title("Raw index test")
	pl.contourf(rm1,phim1,inx1.real,5,cmap=pl.cm.jet)
	pl.contour(rm2,phim2,inx2.real,5,colors='w')

	pl.subplot(212)
	pl.title("Registration test")
	wg.plot(fill=0, sectors=1)
	wg.plot_bitmap(Nx, sectors = 1)

	pl.show()

	Nx=(10,11)

	c1 = wg.get_coord(Nshape=Nx)
	c2 = coordinates.PolarCoord(rrange=(0,wg.rmax+2*c1.dr),
		arange=(-pi/wg.symmetry,pi/wg.symmetry), N=(Nx[0]+2,Nx[1]), border=0)

	#Get raw index and compare to resampled index
	inx1 = wg.index(1.0, coord=c1)
	inx2 = wg.index(1.0, coord=c1, resample=c2)
	print fft.fftshift(inx1, axes=[1]).real
	print inx2.real


if __name__ == "__main__":
	#testing.run_module_suite()
	pass
