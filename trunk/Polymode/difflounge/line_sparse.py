#
# Lined Sparse Class
#
from numpy import *

class Line(object):
	def __init__(self, value=None, offset=None, offsetr=None):
		self.value = asarray(value)
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
		return dot(x.T,self.value).T

	def rmatvec(self, x):
		s = (slice(None),)*self.value.ndim + (newaxis,)*x.ndim
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
		if self.has_key(ii):
			return dict.__getitem__(self,ii)
		else:
			return self.default

class LineSparse(object):
	def __init__(self, shape=(0,0), dtype=complex_):
		self.lines = DefaultDict()
		self.default = None
		self.dtype = dtype
		self.shape = shape

		self.line_range = 0, shape[0]

	def set_size(self,N):
		self.shape = (N,N)
	
	def mirror(self, x, ii=0):
		l = self.lines[ii]
		if l.ndim==1:
			mv = (x.T*l.value).T
		else:
			mv = x*l.value
		return mv

	def matvec(self, x, axis=0):
		"Multiply vector x by implicit matrix"
		mv = zeros(x.shape, self.dtype)
		for ii in range(*self.line_range):
			l = self.lines[ii]
			if l is not None:
				mv[ii] = l.matvec(x[ii-l.offset:ii-l.offset+l.shape[0]])
		return mv

	def rmatvec(self, x):
		"Multiply vector x by adjoint of implicit matrix"
		mv = zeros(x.shape, self.dtype)
		for ii in range(*self.line_range):
			l = self.lines[ii]
			if l is not None:
				xslice = slice(ii-l.offset,ii-l.offset+l.shape[0])
				mv[xslice] += l.rmatvec(x[ii])
		return mv

	def toarray(self):
		A = zeros(self.shape, dtype=self.dtype)
		for ii in range(self.shape[0]):
			l = self.lines[ii]
			if l is not None:
				if l.value.ndim==1:
					A[ii,ii-l.offset:ii-l.offset+l.shape[0]] = l.value
				else:
					A[ii,ii-l.offset:ii-l.offset+l.shape[0]] = l.value[:,0]

		return A


if __name__=="__main__":
	from mathlink import timer
	
	L = DefaultDict()
	L.default = Line([1,-2,1])
	L[0] = Line([-2,1,3,4], 0)
	L[9] = Line([1,-2], 1)

	#Run
	N = 10
	
	AS = LineSparse((N,N))
	AS.lines.default = Line([1,-2,1])
	AS.lines[0] = Line([1,1,-2,1], 0)
	AS.lines[N-1] = Line([1,-2], 1)

	A = AS.toarray()

	x = random.random(N)
	y = AS.matvec(x)
	print "matvec error:", abs(y-dot(A,x)).max()
	
	yr = AS.rmatvec(x)
	print "matvec error:", abs(yr-dot(A.T,x)).max()

	import timer
	
	#Time it
	tick = timer.timer()
	
	tick.start('matvec')
	for ii in range(100):
		y = AS.matvec(x)
	tick.stop('matvec')

	tick.start('rmatvec')
	for ii in range(100):
		y = AS.rmatvec(x)
	tick.stop('rmatvec')

	print tick.report()

	#Try multidimensional arrays
	x = random.random((N,2,5))
	coeff = zeros((4,2,5))
	coeff.T[:] = [1,1,-2,1]
	coeff[:,0,0] = [4,4,-1,1]
	
	N = 10
	AS = LineSparse((N,N))
	AS.lines.default = Line([1,-2,1])
	AS.lines[0] = LineMulti(coeff, 0)
	AS.lines[N-1] = Line([1,-2], 1)

	y2 = AS.matvec(x)
	y2r = AS.rmatvec(x)

	for K in ndindex(x.shape[1:]):
		yd = dot(A, x[(slice(None),)+K])
		print K, abs(yd-y2[(slice(None),)+K]).max()

	for K in ndindex(x.shape[1:]):
		yrd = dot(A.T, x[(slice(None),)+K])
		print K, abs(yrd-y2r[(slice(None),)+K]).max()


