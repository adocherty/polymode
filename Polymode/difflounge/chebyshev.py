#chebsuite.py
# Changelog:
#  6/5/07: Used algorithm by Weideman
#
#
# Todo:
#  
#
#
#

from numpy import size, r_, zeros, empty, eye, diag, pi, tan, arange, real, imag, \
				sin, cos, mod, ndarray, mat, newaxis, sum, absolute
from numpy.fft import fft,ifft,fftshift

class ChebyshevBase:
	def xnodes(self):
		indices = arange(0,self.N)
		x = cos(pi*indices/(self.N-1))*self.L
		return x
	
	def cheb_fft_d(self, corr_x,y,n=1):
		N=len(y)
		theta = r_[0:N-2,0,-N+1:0]

		#fft of y
		yhat = fft(r_[y,flipud(y[1:N-1])])
		yp = -ifft((1j*theta)**n * yhat)[0:N]*corr_x

		#endpoints:
		yp[0] = sum((theta**2*yhat)[0:N])/(N-1)+0.5*(N-1)*yhat[N]
		yp[-1] = -sum(((-1)**theta*theta**2*yhat)[0:N])/(N-1)+0.5*N*(-1)**(N-1)*yhat[N]

		return yp #/self.L

	def diff_1(self, x,y):
		yp = cheb_fft_d(1./sqrt(1-x*x),y,n=1)
		return yp	#/self.L

	def diff_2(self, x,y):
		yp = -cheb_fft_d(1./(1-x*x),y,n=2) \
			+ cheb_fft_d(x/(1-x*x)**1.5,y,n=1)
		return yp	#/self.L


class ChebyshevDense(ChebyshevBase):
	def __init__(self, N, L=1, dtype=float, leftn=0, rightn=0):
		self.N = N			#Number of nodes
		self.L = L			#Rescale to [-L,L]
		self.dtype = dtype
		self.leftn = leftn
		self.rightn = rightn

		self.cache = {}
		
	def diff_matrix(self, dn=1):
		if 'N' not in self.cache or self.N != self.cache['N']:
			self.cache = {}
			self.cache['N'] = self.N
			self.cache['maxdn'] = 0

		if self.cache['maxdn'] < dn:
			self.__construct__(self.N, self.cache['maxdn'], dn)
			self.cache['maxdn'] = dn
		return self.cache[dn]/self.L

	def __construct__(self, N, dn_from, dn_to):
		indices = arange(0,N)
		theta = pi*indices/(N-1)
		
		#xi-xj calculated with trig identity
		Dx = 2*sin(0.5*(theta+theta[:,newaxis])) * sin(0.5*(theta-theta[:,newaxis]))
		Dx[indices, indices] = 1
		IDx = 1/Dx
		IDx[indices, indices] = 0
		
		#Off-diagonal coefficient cij = (-1)**(i+j)
		#i=0 or i=N-1 gives cij = 2*(-1)**(i+j)
		c = (-1.)**indices; c[0]*=2; c[-1]*=2

		#Storage
		if dn_from==0:
			self.cache[0] = eye(N)
		
		for ii in range(dn_from, dn_to):
			D = self.cache[ii]

			#Off-diagonal matrix
			Dnew = (ii+1)*(c[:,newaxis]*D.diagonal()[:,newaxis]/c - D)*IDx

			#Correct diagonal
			Dnew[indices, indices] = -sum(Dnew, axis=1)
			
			Dnew[absolute(Dnew)<1e-12] = 0
			
			self.cache[ii+1] = Dnew

if __name__ == "__main__":
	import timer

	## Test and time!
	tick = timer.timer()

	N=500
	cd=ChebyshevDense(N)

	tick.start("New Method")
	cd.diff_matrix(5)
	tick.stop()

	tick.print_list()
