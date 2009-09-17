#ifndef _AMOS_H
#define _AMOS_H

#include <iostream>
#include <complex>
#include <cmath>
#include "cephes/mconf.h"

#define FNAME(x) x##_

/*
int cairy_(double *, int *, int *, doublecomplex *, int *, int *);
int cbiry_(doublecomplex *, int *, int *, doublecomplex *, int *, int *);
int cbesi_(doublecomplex *, double *, int *, int *, doublecomplex *, int *, int *);
int cbesj_(doublecomplex *, double *, int *, int *, doublecomplex *, int *, int *);
int cbesk_(doublecomplex *, double *, int *, int *, doublecomplex *, int *, int *);
int cbesy_(doublecomplex *, double *, int *, int *, doublecomplex *, int *, doublecomplex *, int *);
int cbesh_(doublecomplex *, double *, int *, int *, int *, doublecomplex *, int *, int *);
*/

/*
External fortran named functions:

Hankel function
ZBESH(ZR, ZI, FNU, KODE, M, N, CYR, CYI, NZ, IERR)
ZR,ZI	- Z=CMPLX(ZR,ZI), Z.NE.CMPLX(0.0D0,0.0D0),
		-PT.LT.ARG(Z).LE.PI
FNU    - ORDER OF INITIAL H FUNCTION, FNU.GE.0.0D0
KODE   - A PARAMETER TO INDICATE THE SCALING OPTION
	KODE= 1  RETURNS
		     CY(J)=H(M,FNU+J-1,Z),   J=1,...,N
		= 2  RETURNS
		     CY(J)=H(M,FNU+J-1,Z)*EXP(-I*Z*(3-2M))
		          J=1,...,N  ,  I**2=-1
M      - KIND OF HANKEL FUNCTION, M=1 OR 2
N      - NUMBER OF MEMBERS IN THE SEQUENCE, N.GE.1

CY(I)=EXP(-MM*Z*I)*H(M,FNU+J-1,Z)       MM=3-2*M,   I**2=-1.
*/
namespace amos
{
	using namespace std;
	
	typedef complex<double> complex_t;
	
	extern "C"
	{
	int FNAME(zairy)(double*, double*, int*, int*, double*, double*, int*, int*);
	int FNAME(zbiry)(double*, double*, int*, int*, double*, double*, int*, int*);
	int FNAME(zbesi)(double*, double*, double*, int*, int*, double*, double*, int*, int*);
	int FNAME(zbesj)(double*, double*, double*, int*, int*, double*, double*, int*, int*);
	int FNAME(zbesk)(double*, double*, double*, int*, int*, double*, double*, int*, int*);
	int FNAME(zbesy)(double*, double*, double*, int*, int*, double*, double*, int*, double*, double*, int*);
	int FNAME(zbesh)(double*, double*, double*, int*, int*, int*, double*, double*, int*, int*);
	}

	//typedef complex<double> complex<double> ;

	complex<double> rotate(complex<double> z, double v)
	{
		double c = cos(v * M_PI);
		double s = sin(v * M_PI);
		complex<double> w(real(z)*c - imag(z)*s, real(z)*s + imag(z)*c);
		return w;
	};

	complex<double> besselJ(double v, complex<double> z)
	{
	  int n = 1;
	  int kode = 1;
	  int nz, ierr;
	  int sign = 1;

	  double z_re = z.real();
	  double z_im = z.imag();
	  
	  double ans_re, ans_im;
	  
	  if (v < 0)
	  {
		v = -v;
		sign = -1;
	  }
	  
	  FNAME(zbesj)(&z_re, &z_im, &v,  &kode, &n, &ans_re, &ans_im, &nz, &ierr);
	  
	  complex<double> ans(ans_re, ans_im);
	  
	  if (sign == -1)
		ans = rotate(ans, v);
	  return ans;
	};

	template<typename T>
	complex<double> besselJp(int m, double v, T z)
	{
		return besselJ(m,v-1,z) - (v/static_cast<complex<double> >(z))*besselJ(m,v,z);
	};

	complex<double> besselK(double v, complex<double> z)
	{
	  int n = 1;
	  int kode = 1;
	  int nz, ierr;

	  double z_re = z.real();
	  double z_im = z.imag();
	  
	  double ans_re, ans_im;
	  
	  if (v < 0) v = -v;
	 
	  FNAME(zbesk)(&z_re, &z_im, &v,  &kode, &n, &ans_re, &ans_im, &nz, &ierr);

	  return complex<double>(ans_re, ans_im);
	};

	// Hankel functions
	complex<double> hankel(int m, double v, complex<double> z)
	{
	  int n = 1;
	  int kode = 1;		//If this is 2 the result is exponentially scaled
	  int sign = 1;
	  int nz, ierr;

	  double z_re = z.real();
	  double z_im = z.imag();
	  
	  double ans_re, ans_im;
	  
	  // zbesh only accepts positive v
	  if (v < 0)
	  { v = -v; sign = -1; }

	  //Access the AMOS bessel routine	 
	  FNAME(zbesh)(&z_re, &z_im, &v, &kode, &m, &n, &ans_re, &ans_im, &nz, &ierr);
	  complex<double> ans(ans_re, ans_im);
	  
	  //Use relation H^(1)_-m(x) = exp(m pi i) H^(1)_m(x)
	  //        and  H^(2)_-m(x) = exp(-m pi i) H^(2)_m(x)
	  if (sign==-1)
	  	if (m==1)      ans = rotate(ans, v);
		else if (m==2) ans = rotate(ans, -v);

	  return ans;
	};

	template<typename T>
	complex<double> hankelp(int m, double v, T z)
	{
		return hankel(m,v-1,z) - (v/static_cast<complex<double> >(z))*hankel(m,v,z);
	};

	//Shortcuts for the hankel function
	template<typename T>
	complex<double> hankel1(double v, T z)
	{
		return hankel(1, v, z);
	};
	template<typename T>
	complex<double> hankel2(double v, T z)
	{
		return hankel(2, v, z);
	};
	template<typename T>
	complex<double> hankel1p(double v, T z)
	{
		return hankel1(v-1,z) - hankel1(v,z)*v/z;
	};
	template<typename T>
	complex<double> hankel2p(double v, T z)
	{
		return hankel2(v-1,z) - hankel2(v,z)*v/z;
	};



}

#endif

