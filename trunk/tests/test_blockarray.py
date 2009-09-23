#Unit test for VWE_Action

from __future__ import division

from numpy import *
from numpy.testing import *

#--------------------------------------------------------------------------
import sys,os
sys.path.append(os.pardir)
#--------------------------------------------------------------------------

from Polymode.mathlink.blockarray import *

class test_blockarray(TestCase):
	blocktest = array([[ 0.,  0.,  1.,  2.,  3.,  4.],
		 [ 0.,  0.,  5.,  6.,  7.,  8.],
		 [ 1.1,  2.1,  3.1,  4.1,  5.1, 6.1],
		 [ 7.1,  8.1,  9.1,  10.1,  11.1,  12.1],
		 [ -1.,  -2.,  -3.,  -4.,  0., 0.],
		 [ -5.,  -6.,  -7.,  -8.,  0.,  0.]])
	densetest = array([[  1.,    2.,    3.,    4.,    0.,    0. ],
			 [  5.,    6.,    7.,    8.,    0.,   0. ],
			 [  1.1,  2.1,   3.1,   4.1,   5.1,   6.1],
			 [  7.1,  8.1,   9.1,  10.1,  11.1,  12.1],
			 [  0.,    0.,   -1.,   -2.,   -3.,   -4.],
			 [  0.,    0.,   -5.,   -6.,   -7.,   -8.]])
	
	def test_dense(self):
		self.__check_dense__(int)
		self.__check_dense__(float)
		self.__check_dense__(complex)

	def test_matvec(self):
		self.__check_matvec__((20,3),(5,5),complex)
		self.__check_matvec__((20,3),(5,5),float)
		self.__check_matvec__((20,5),(4,4),complex)
		self.__check_matvec__((20,7),(8,8),complex)

	def test_misc(self):
		self.__check_misc__((20,3),(8,8),complex)
		self.__check_misc__((20,3),(8,8),float)
		
		#Currently doesn't work for more than 3 columns...
		#self.__check_misc__((4,5),(2,2),int)


	def __check_dense__(self, dtype):
		print "Dense: Checking %s" % dtype
		mshape=(3,3); blockshape=(2,2)
		ba = BlockArray(mshape, blockshape, dtype)
		ba[:] = self.blocktest
		dense = self.densetest.astype(dtype)

		assert_almost_equal(ba.toarray(),dense)
		
	def __check_matvec__(self, mshape, blockshape, dtype):
		print "Matvec: Checking size:%s, block:%s for %s" % (mshape,blockshape,dtype)
		ba = BlockArray(mshape, blockshape, dtype)
		ba[:] = 100*random.random(ba.shape)
		
		#Check left and right vector multiplies
		x=random.rand(ba.shape[0])
		dense = ba.toarray()
		assert_almost_equal(ba.matvec(x), dot(dense,x))
		assert_almost_equal(ba.rmatvec(x), dot(dense.T,x))

	def __check_misc__(self, mshape, blockshape, dtype):
		print "Misc: Checking size:%s, block:%s for %s" % (mshape,blockshape,dtype)
		ba = BlockArray(mshape, blockshape, dtype)
		ba[:] = 100*random.random(ba.shape)
		
		#Check block transpose
		assert_almost_equal(ba.transpose().toarray(), ba.toarray().transpose())



class test_blocklu(TestCase):
	'''
	Check BlockLU Functionality
	'''
	def check_blocklu(self):
		#Only check tridiagonal
		self.__check__((100,3),(5,5),float)
		self.__check__((100,3),(5,5),complex)

	def __check__(self, shape, blockshape, dtype):
		print "BlockLU: Checking size:%s, block:%s for %s" % (shape,blockshape,dtype)
		ba = BlockArray(shape, blockshape, dtype)
		ba[:] = 100*random.random(ba.shape)
		
		#Create LU Decomposition
		lutri = TriBlockLU( overwrite=False )
		lutri.lu(ba)
		
		x = random.rand(ba.shape[0]).astype(dtype)

		#Check standard solve:
		y = lutri.solve(x)
		assert_almost_equal(ba.matvec(y), x)

		#Check transpose solve:
		yt = lutri.solve_transpose(x)
		assert_almost_equal(ba.rmatvec(yt), x)

		if hash(ba)==hash(lutri.Alu):
			print "Warning: lutri is overwriting A"

		#Check updates
		lu = lutri.lu_update(ba, uprows=3)
		yup = lutri.solve(x)
		assert_almost_equal(ba.matvec(yup), x)

		ytup = lutri.solve_transpose(x)
		assert_almost_equal(ba.rmatvec(ytup), x)
		
if __name__ == '__main__':
    run_module_suite()
