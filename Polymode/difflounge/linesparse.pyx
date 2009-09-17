from __future__ import division
import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np

#The python complex double
cdef extern from "complexobject.h":
	ctypedef struct Py_complex:
		double real
		double imag

	ctypedef class __builtin__.complex [object PyComplexObject]:
		cdef Py_complex cval

#Wrap the c++ std::complex<double> class
cdef extern from "complex":
	ctypedef struct ccomplex "std::complex<double>":
		double real()
		double imag()
	ccomplex *new_ccomplex "new std::complex<double>" (double x, double y)
	void del_complex "delete" (ccomplex *c)

#def inc1_cfloat(np.ndarray[np.cfloat_t] arr):
#.   arr[1].real += 1
#    arr[1].imag += 1

# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.cdouble

# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.cdouble_t DTYPE_t

# The builtin min and max functions works with Python objects, and are
# so very slow. So we create our own.
#  - "cdef" declares a function which has much less overhead than a normal
#    def function (but it is not Python-callable)
#  - "inline" is passed on to the C compiler which may inline the functions
#  - The C type "int" is chosen as return type and argument types
#  - Cython allows some newer Python constructs like "a if x else b", but
#    the resulting C file compiles with Python 2.3 through to Python 3.0 beta.
cdef inline int int_max(int a, int b): return a if a >= b else b
cdef inline int int_min(int a, int b): return a if a <= b else b

#Wrap indices around, i'm sure there's a better way!
cdef inline int c_mod(int offset,int size):
	cdef int n = offset%size
	return n if n>=0 else n+size

cdef inline np.cdouble_t cc_add(np.cdouble_t a, np.cdouble_t b):
	return np.cdouble_t(a.real+b.real, a.imag+b.imag)

cdef inline np.cdouble_t cc_mul(np.cdouble_t a, np.cdouble_t b):
	return np.cdouble_t(a.real*b.real-a.imag*b.imag, a.real*b.imag+a.imag*b.real)


# "def" can type its arguments but not have a return type. The type of the
# arguments for a "def" function is checked at run-time when entering the
# function.
#
# The arrays f, g and h is typed as "np.ndarray" instances. The only effect
# this has is to a) insert checks that the function arguments really are
# NumPy arrays, and b) make some attribute access like f.shape[0] much
# more efficient. (In this example this doesn't matter though.)

def c_dot(np.ndarray[np.cdouble_t] x, np.ndarray[np.cdouble_t, , ndim=2] y):
	#Check type
	assert x.dtype == np.cdouble and y.dtype == np.cdouble

	#Size of array
	cdef int xsize = x.shape[0]
	cdef int ysize = y.shape[0]

	cdef np.cdouble_t d = np.cdouble_t(0,0)
	#cdef np.ndarray[np.cdouble_t] h = np.zeros(xsize, dtype=np.cdouble)

	cdef int xi
	for xi in range(xsize):
		d = cc_add(d, cc_mul(x[xi], y[xi]))
	return complex(d.real, d.imag)

def offset_dot(np.ndarray[np.cdouble_t] x, np.ndarray[np.cdouble_t] y, int offset=0):
	#Check type
	assert x.dtype == np.cdouble and y.dtype == np.cdouble

	#Size of array
	cdef int xsize = x.shape[0]
	cdef int ysize = y.shape[0]

	#ysize can be < xsize
	assert ysize<=xsize

	cdef np.cdouble_t d = np.cdouble_t(0,0)
	#cdef np.ndarray[np.cdouble_t] h = np.zeros(xsize, dtype=np.cdouble)

	cdef int yi
	for yi in range(ysize):
		d = cc_add(d, cc_mul(x[c_mod(yi + offset,xsize)], y[yi]) )
	return complex(d.real, d.imag)

class Line:
	def __init__(self, value, offset=None, offsetr=None):
		self.value = np.asarray(value)
		self.shape = self.value.shape
		self.ndim = self.value.ndim
		
		if offset is None:
			offset = self.shape[0]//2
		if offsetr is not None:
			offset = self.shape[0]-offsetr
		self.offset = offset

	def __str__(self):
		return "<Line: %s @ %i>" % (self.value, self.offset)

	def matvec(self, x):
		return np.dot(x.T,self.value).T

	def rmatvec(self, x):
		s = (slice(None),)*self.value.ndim + (np.newaxis,)*x.ndim
		return self.value[s]*x

class LineMulti(Line):
	"Line object for multiple dimensions"
	def matvec(self, x):
		return sum(x*self.value, axis=0)
	
	def rmatvec(self, x):
		return self.value*x

class DefaultDict(dict):
	"Dictionary with default (fallback) value"
	def __init__(self):
		self.default = None
		dict.__init__(self)
		
	def __getitem__(self, ii):
		try:
			return dict.__getitem__(self,ii)
		except KeyError:
			return self.default

class LineSparse(object):
	def __init__(self, shape=(0,0)):
		self.lines = DefaultDict()
		self.default = None
		self.shape = shape

		self.line_range = 0, shape[0]

	def set_size(self, N):
		self.shape = (N,N)

	def mirror(self, x, ii=0):
		l = self.lines[ii]
		if l.ndim==1:
			mv = (x.T*l.value).T
		else:
			mv = x*l.value
		return mv

	def matvec(self, np.ndarray[DTYPE_t, ndim=2] x, int axis=0):
		"Multiply vector x by implicit matrix"
		cdef np.ndarray[DTYPE_t, ndim=2] mv = np.zeros([x.shape[0], x.shape[1]], dtype=DTYPE)
		cdef object l
		cdef int offset
		cdef int ii
		
		for ii in range(self.line_range[0], self.line_range[1]):
			l = self.lines[ii]
			if l is not None:
				offset = ii-l.offset
				#mv[ii] = offset_dot(x, l.value, offset)
				mv[ii] = np.dot(x[offset:offset+l.shape[0]].T,l.value).T
		return mv

	def rmatvec(self, x):
		"Multiply vector x by adjoint of implicit matrix"
		mv = np.zeros(x.shape, DTYPE)
		for ii in range(*self.line_range):
			l = self.lines[ii]
			if l is not None:
				xslice = slice(ii-l.offset,ii-l.offset+l.shape[0])
				mv[xslice] += l.rmatvec(x[ii])
		return mv

	def toarray(self):
		A = np.zeros(self.shape, DTYPE)
		for ii in range(self.shape[0]):
			l = self.lines[ii]
			if l is not None:
				if l.value.ndim==1:
					A[ii,ii-l.offset:ii-l.offset+l.shape[0]] = l.value
				else:
					A[ii,ii-l.offset:ii-l.offset+l.shape[0]] = l.value[:,0]

		return A

