

from numpy import size, r_, zeros, pi, tan, arange, real, imag, empty, sin, mod, ndarray, eye
from numpy.fft import fft,ifft,fftshift

# Generate Fourier differentiation matricies
# Two routines to do so: standard and fourier methods
# 
# 7/1/07 - Fixed now routines work for all N & power
#
# Todo:
# * Clean up methods .. should be one. Do they work?
# * Special methods for powers 1 & 2
# * Save last k DMs and lookup N & power

## Matlab-like Toeplitz function
def toeplitz(D, col_left, col_right = 0):
	if not issubclass(type(col_right), np.ndarray):
		col_right = col_left
	
	N = len(col_left)
	if N!=len(col_right):
		raise IndexError("Lengths of col_left & col_right not the same")
	if D.shape!=(N,N):
		raise IndexError("Matrix not correct shape")

	for j in np.arange(0,N):
		D[j] = np.append(col_left[j::-1], col_right[1:N-j])

def fourier_space(N, symmetric = False):
		if symmetric and N%2==0:
			ks = r_[0:N//2,0,-N//2+1:0]
		else:
			ks = r_[0:(N+1)//2,-(N//2):0]
		return ks
			
class FourierAction:
	def xnodes(self):
		pass

	def diff_n(self, y, dn=1, axis=0):
		assert isinstance(y, ndarray), "Not a NumPy array type"
		
		k = fourier_space(y.shape[axis])
		yp = empty(y.shape, dtype = y.dtype)
	
		# Differentition in Fourier space, selecting correct axis
		yp.swapaxes(axis, y.ndim)[:] = ifft((1j*k)**dn * fft(y.swapaxes(axis, y.ndim), axis=y.ndim), axis=y.ndim)
		return yp
	
	def __call__(self, y, dn=1, axis=0):
		return self.diff_n(y,dn,axis)

class FourierDense:
	def __init__(self, N, dtype=float):
		self.N = N
		self.dtype = dtype

	def dft_matrix(self):
		F = empty((self.N,self.N),dtype=self.dtype)
		I = eye(self.N)
		for i in range(self.N):
			F[i] = fft(I[i])
		return F

	def idft_matrix(self):
		IF = empty((self.N,self.N),dtype=self.dtype)
		I = eye(self.N)
		for i in range(self.N):
			IF[i] = ifft(I[i])
		return IF
	
	def diff_matrix(self, power=1):
		indices = arange(0,self.N)
		
		#Choose symmetric k
		if (mod(self.N,2)==0) and (mod(power,2)==1):
			k = keven = r_[0:self.N//2,0,-self.N//2+1:0]
		else:
			k = kodd = r_[0:(self.N+1)//2,-(self.N//2):0]

		kdelta = zeros(self.N, dtype=self.dtype)
		kdelta[0] = 1
		col=real(ifft( (1j*k)**power * fft(kdelta) ))

		if mod(power,2)==0:
			colx = col
		else:
			col[0]=0
			colx = -col
	
		D = empty((self.N,self.N),dtype=self.dtype)
		toeplitz(D,col,colx)
		return D
	
