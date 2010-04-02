#ifndef _BESSELRATIOS_H
#define _BESSELRATIOS_H

#include <iostream>
#include <complex>
#include <cmath>

#include "amos.h"

#define FNAME(x) x##_

namespace besselratio
{
	using namespace std;
	
	typedef complex<double> complex_t;

	double RATIO_MMAX = 10;
	int MAXITER = 100;
	
	//Series expansion for the bessel ration J_{m+1}(x)/J_{m-1}(x)
	template<typename T>
	T  besselj_ratio_series(double m, T z)
	{
		T S;
		double sign = 1;
		double v = abs(m);
		T Sup = 0;
		T Slow = 0;
		T Gns = 1;
		T Gs = 1;
		
		double tol = 1e-10;
		double s = 0;
		T correction = 1;
		while( (abs(correction)>tol) && (s<MAXITER) )
		{
			Gns *= (v+s-1.0);
			S = sign * pow(0.5*z,2*s)/(Gs*Gns);
			Gs *= (s+1.0);
			sign *= -1.0;
			
			Sup += S/(v+s)/(v+s+1.0);
			Slow += S;
			s += 1;
			
			correction = S;
		}

		return pow(0.5*z, 2) * Sup/Slow;
	};

	//Continues fraction expression  for the bessel ration J_{m+1}(x)/J_{m-1}(x)
	//From Lentz (1979)
	template<typename T>
	T  besselj_ratio_cnf(double m, T z)
	{
		T fn, C0, C1;
		double v = abs(m)+1;
		double s = 1;
		C1 = 0;
		C0 = 1.0/(2*v/z);
		fn = C0;

		double tol = 1e-12;
		T correction = 0;
		int niter = 0;
		while( abs(correction-1.0)>tol && (niter<MAXITER) )
		{
			v += 1;
			s *= -1;
			C0 = 1.0/(s*2*v/z + C0);
			C1 = 1.0/(s*2*v/z + C1);
			
			correction = C0/C1;
			fn *= correction;
			niter += 1;
		}

		return fn/(2*m/z-fn);
	};

	//Hankel function ratio  H_{m-1}(x)/H_{m+1}(x)
	//From recurrance relations
	template<typename T>
	T  hankel1_ratio_recurrance(double m, T z)
	{
		double v = abs(m);
		T Hfn = amos::hankel1(0.0,z)/amos::hankel1(1.0,z);

		double tol = 1e-12;
		T correction = 0;
		for(double s=0; s<v-1; s++)
		{
			//Recurrance relation
			Hfn = 1.0/(2.0*(s+1)/z - Hfn);
		}

		return Hfn/(2*m/z-Hfn);
	};

	template<typename T>
	T  hankel1_ratio(double m, T z)
	{
		double ma = abs(m);
		T hmhp, HR;

		//For small m amos functions are fine, only use ratio routines for large m		
		if(ma<RATIO_MMAX)
			HR = ma - z*amos::hankel1(ma+1,z)/amos::hankel1(ma,z);
		else
		{
			hmhp = hankel1_ratio_recurrance(ma, z);
			HR = -ma*(1.0 - hmhp)/(1.0 + hmhp);
		}
		return HR;
	};

	//Calculate x J'_m(x)/J_m(x)
	template<typename T>
	T  besselj_ratio(double m, T z)
	{
		double ma = abs(m);
		T jpjm, BR;
		
		if(ma<RATIO_MMAX)
			BR = ma - z*amos::besselJ(ma+1,z)/amos::besselJ(ma,z);
		else
		{
			jpjm = besselj_ratio_cnf(ma, z);
			BR = ma*(1.0 - jpjm)/(1.0 + jpjm);
		}
		return BR;
	};

}

#endif

