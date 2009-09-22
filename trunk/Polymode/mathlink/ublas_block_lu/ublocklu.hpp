#ifndef UBLASBLOCKLU_H
#define UBLASBLOCKLU_H

#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/operation.hpp>

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

//Boost bindings to atlas
#include <boost/numeric/bindings/atlas/cblas2.hpp>
#include <boost/numeric/bindings/atlas/cblas3.hpp>
#include <boost/numeric/bindings/atlas/clapack.hpp>
#include <boost/numeric/bindings/traits/ublas_matrix.hpp>
#include <boost/numeric/bindings/traits/std_vector.hpp>

//Boost python:
#include <Python.h>
#include <boost/python.hpp>
#include <boost/python/numeric.hpp>

//Numpy
#include <numpy/noprefix.h>
#include <numpy/arrayobject.h>

#include "numimport.hpp"

const char UBLASBLOCKLU_VERSION[] = "v0.1.0";

namespace numpy = boost::python::numeric;
namespace ublas = boost::numeric::ublas;
namespace atlas = boost::numeric::bindings::atlas;

template<class T>
void print_matrix(const ublas::matrix_expression<T>& Ae)
{
	T A = Ae();
	
	std::cout << "[";
    for(unsigned i = 0; i < A.size1(); ++i) {
	    for(unsigned j = 0; j < A.size2(); ++j)
	    	printf("%8.4g ", A(i,j));
	    if( i<(A.size1()-1) ) std::cout << "\n ";
	}
	std::cout << " ]\n";
};

typedef std::complex<double> complex_t;

/*Wrapper for block matrix
  To wrap existing data with ublas matricies choose one of these two:
  neither seem to copy the data
  
	vector<double, array_adaptor<double> > v(5);
	v.data().resize(5, data);

    shallow_array_adaptor<double> shared(5, data);
    vector<double, shallow_array_adaptor<double> > shared_vec(5, shared);
*/
template<class T>
class block_array_wrapped
{
private:
	typedef typename ublas_extern<T>::adaptor_type ad_t;
	typedef typename ublas_extern<T>::external_matrix_type exmat_t;

	ad_t* shared_adaptor_ptr;
	exmat_t* wrapped_array_ptr;
	unsigned nrows;
	unsigned ncols;
	unsigned nblock;
	
public:
	//Contructors
	block_array_wrapped() : shared_adaptor_ptr(0), wrapped_array_ptr(0) {}

	block_array_wrapped(unsigned nrows_in, unsigned ncols_in, unsigned nblock_in, T* data)
	{ 
		nblock=nblock_in;		nrows=nrows_in;		ncols=ncols_in;
		shared_adaptor_ptr = new ad_t(nrows_in*ncols_in*nblock_in*nblock_in, data);
		wrapped_array_ptr = new exmat_t(nrows*nblock, ncols*nblock, *shared_adaptor_ptr);
	}

	block_array_wrapped( block_array_wrapped<T>& M )
	{
		nrows = M.nrows; ncols = M.ncols; nblock = M.nblock;
		shared_adaptor_ptr = new ad_t( *M.shared_adaptor_ptr );
		wrapped_array_ptr = new exmat_t(nrows*nblock, ncols*nblock, *shared_adaptor_ptr);
	}

	~block_array_wrapped()
	{
		delete_arrays();
	}

	void delete_arrays()
	{
		if( shared_adaptor_ptr>0 ) delete shared_adaptor_ptr;
		if( wrapped_array_ptr>0 ) delete wrapped_array_ptr;
	}

	void assign( block_array_wrapped<T>& M )
	{
		delete_arrays();
		nrows = M.nrows; ncols = M.ncols; nblock = M.nblock;
		shared_adaptor_ptr = new ad_t( *M.shared_adaptor_ptr );
		wrapped_array_ptr = new exmat_t(nrows*nblock, ncols*nblock, *shared_adaptor_ptr);
	}

	void assign(unsigned nrows_in, unsigned ncols_in, unsigned nblock_in, T* data)
	{ 
		nblock=nblock_in;		nrows=nrows_in;		ncols=ncols_in;
		shared_adaptor_ptr = new ad_t(nrows_in*ncols_in*nblock_in*nblock_in, data);
		wrapped_array_ptr = new exmat_t(nrows*nblock, ncols*nblock, *shared_adaptor_ptr);
	}

	exmat_t& matrix()
	{ return *wrapped_array_ptr; }
	
	const unsigned get_nrows() { return nrows; }
	const unsigned get_ncols() { return ncols; }
	const unsigned get_nblock() { return nblock; }
	
	ublas::vector<unsigned> shape()
	{
		ublas::vector<unsigned> shape(2);
		shape[0] = matrix().size1();
		shape[1] = matrix().size2();
		return shape;
	}
	
	template<class AE>
	void set_block(unsigned i, unsigned j, const ublas::matrix_expression<AE>& A)
	{
		ublas::subrange(matrix(), i*nblock, (i+1)*nblock, j*nblock, (j+1)*nblock) = A;
	}
	ublas::matrix_range<exmat_t> get_block(unsigned i, unsigned j)
	{
		return ublas::subrange(matrix(), i*nblock, (i+1)*nblock, j*nblock, (j+1)*nblock);
	}

};

/****************************************************************************************

Class to perform actual LU

****************************************************************************************/
template<class T>
class tridiagonal_block_lu
{
private:
	bool overwrite;
	T shift;
	ublas::matrix<std::size_t> pivots_m;
	block_array_wrapped<T> LU;

public:
	tridiagonal_block_lu(const tridiagonal_block_lu<T>& lu ) {};
	tridiagonal_block_lu(bool overwrite_in = true )
		{ overwrite = overwrite_in; }
	bool isoverwrite() { return overwrite; }
	
	//copy the pivot information for row irow into a new pivots matrix
	std::vector<int> get_pivots(unsigned irow)
	{
		std::vector<int> pivots(pivots_m.size2());
		for(unsigned i=0; i<pivots.size(); ++i)
			pivots[i] = pivots_m(irow,i);
		return pivots;
	}

	void store_pivots(std::vector<int> pivots, unsigned irow)
	{
		for(unsigned i=0; i<pivots.size(); ++i)
			pivots_m(irow,i) = pivots[i];
	}

	//Perform Block LU
	void lu(block_array_wrapped<T>& A, T shift);
	void lu_update(block_array_wrapped<T>& A, int nrowsup);

	//Back subsistution
	template<class M, class S>
	void solve(ublas::matrix<T,M,S>& x);

	//Back subsistution of transpose problem
	template<class M, class S>
	void solve_transpose(ublas::matrix<T,M,S>& x);

	//Python interface
	void py_lu(numpy::array& Ain, unsigned nblock, T shift);
	void py_lu_update(numpy::array& Aupdate, int uprows);
	void py_solve(numpy::array& x);
	void py_solve_transpose(numpy::array& x);
};

template<class T>
void tridiagonal_block_lu<T>::lu(block_array_wrapped<T>& A, T shift_in=0)
{
    using namespace ublas;
	using namespace std;

	unsigned nrows = A.get_nrows();
	unsigned ncols = A.get_ncols();
	unsigned nblock = A.get_nblock();

	shift = shift_in;

	//Check size
	if( ncols!=3 ) cout << "Error: Only tridiagonal block matricies are supported!\n";

	pivots_m.resize(nrows, nblock);			//Resize pivots storage

	//Overwrite current matrix?
	//if( overwrite )
	
	//Temporary storage of blocks
	matrix<T, column_major> a,b,c;
	matrix<T, column_major> Ishift = shift * identity_matrix<T>(nblock);
	
	//Shift first block
	b = A.get_block(0, 1) - Ishift;
	A.set_block(0, 1, b);

	//Loop over block rows
	for( unsigned i=0; i<nrows-1; ++i)
	{
		a = A.get_block(i+1, 0);
		c = A.get_block(i, 2);
		
		//real permulation matrix for LU routines
		//Needs to be here to reset status of pivots each time (or for loop)
		std::vector<int> pivots(nblock);

		//LU factorize on block
		atlas::lu_factor(b, pivots);
		
		int isnan_flag=0;
		for(unsigned xx=0; xx<nblock; xx++ )
			for(unsigned yy=0; yy<nblock; yy++ )
				if( isnan(real(static_cast<complex_t>(b(xx,yy)))) )
					isnan_flag = 1;
		if( isnan_flag )
		{
			cout << "Nan detected in LU, row=" << i << endl;
		}
			
		//Solve anew = a b^-1
		a = trans(a);
		atlas::lu_substitute_transpose(b, pivots, a);
		a = trans(a);
		
		store_pivots(pivots, i);			//Store pivots
		
		//Store changed blocks
		A.set_block(i,1,b);					//blu
		A.set_block(i+1,0,a);				//anew

		//bnew
		b = A.get_block(i+1, 1);
		noalias(b) -= Ishift;
		atlas::gemm(-1, a, c, 1, b);
		//qnoalias(b) -= prod(a, c);

		//Copy over block c, if not overwriting
		//if( !self.overwrite ) LUinternal.set_block(i,2,c);
	}

	//Final LU decomp of last block
	std::vector<int> pivots(nblock);
	atlas::lu_factor(b, pivots);
	store_pivots(pivots, nrows-1);
	A.set_block(nrows-1, 1, b);
};

template<class T>
void tridiagonal_block_lu<T>::lu_update(block_array_wrapped<T>& Aupdate, int uprows)
{
    using namespace ublas;
	using namespace std;

	unsigned nrowsup = Aupdate.get_nrows();
	unsigned nrows = LU.get_nrows();
	unsigned nblock = LU.get_nblock();
	
	//Temporary storage of blocks
	matrix<T, column_major> a,b,c;
	matrix<T, column_major> Ishift = shift * identity_matrix<T>(nblock);

	std::vector<int> pivots(nblock);
	
	//b[i-1] is already LU decomposed from original decomposition
	b = LU.get_block(nrows-uprows-1, 1);
	pivots = get_pivots(nrows-uprows-1);
	for(int i=-uprows; i<0; ++i)
	{
		a = Aupdate.get_block(nrowsup+i,0);
		if( i==-uprows )
			c = LU.get_block(nrows+i-1,2);			//Get c[i-1] from original matrix at start row
		else
			c = Aupdate.get_block(nrowsup+i-1,2);	//Then from update at sucessive rows
		
		//anew[i] bnew[i-1] = a[i]
		a = trans(a);
		atlas::lu_substitute_transpose(b, pivots, a);
		a = trans(a);

		//bnew[i] =  b[i] - a[i-1]c[i]
		b = Aupdate.get_block(nrowsup+i, 1);
		noalias(b) -= Ishift;
		atlas::gemm(-1, a, c, 1, b);

		atlas::lu_factor(b, pivots);

		//Store changed blocks
		store_pivots(pivots, nrows+i);
		LU.set_block(nrows+i,0,a);					//anew
		LU.set_block(nrows+i,1,b);					//blu
		LU.set_block(nrows+i-1,2,c);
	}

};


/*
	Solve LU y = x
	with x.shape = (nrows, nblock)
*/
template<class T>
template<class M, class S>
void tridiagonal_block_lu<T>::solve(ublas::matrix<T,M,S>& x)
{
    using namespace ublas;
	using namespace std;

	unsigned nrows = LU.get_nrows();
	unsigned nblock = LU.get_nblock();

	//Temporaries
	matrix<T, column_major> b(nblock,nblock);
	matrix<T, column_major> xrm(nblock,1);
	ublas::vector<T> xr(nblock);

	//Check matrix size
	if( (x.size1()!=nrows) || (x.size2()!=nblock) )
		cout << "Error: x shape is not correct" << endl;

	//Forward substitution pass
	// dnew[0] = d[0]
	// dnew[i] = d[i]-anew[i].dnew[i-1]
	for(unsigned i=0; i<nrows-1; ++i)
	{
		xr = row(x,i+1);
		atlas::gemv(-1, LU.get_block(i+1,0), row(x,i), 1, xr);
		row(x,i+1) = xr;
		//row(x,i+1) -= prod(LU.get_block(i+1,0), row(x,i));
	}
	
	//Backward substitution
	for(int i=nrows-1; i>=0; --i)
	{
		if(i<int(nrows-1))
		{	xr = row(x,i);
			atlas::gemv(-1, LU.get_block(i,2), row(x,i+1), 1, xr);
			row(x,i) = xr;
			//row(x,i) -= prod(LU.get_block(i,2), row(x,i+1));
		}
		b = LU.get_block(i,1);

		column(xrm,0) = row(x,i);
		atlas::lu_substitute(b, get_pivots(i), xrm);
		row(x,i) = column(xrm,0);
	}

};

template<class T>
template<class M, class S>
void tridiagonal_block_lu<T>::solve_transpose(ublas::matrix<T,M,S>& x)
{
    using namespace ublas;
	using namespace std;

	unsigned nrows = LU.get_nrows();
	unsigned nblock = LU.get_nblock();

	//Temporaries
	matrix<T, column_major> b(nblock,nblock);
	matrix<T, column_major> xrm(nblock,1);
	ublas::vector<T> xr(nblock);

	//Check matrix size
	if( (x.size1()!=nrows) || (x.size2()!=nblock) )
		cout << "Error: x shape is not correct" << endl;

	//Backward substitution
	// b[0].T dnew[0] = d[0]
	// b[i].T dnew[i] = d[i] - c[i-1].T dnew[i-1]
	for(unsigned i=0; i<nrows; ++i)
	{
		if( i>0 )
		{	xr = row(x,i);
			atlas::gemv(CblasTrans, -1, LU.get_block(i-1,2), row(x,i-1), 1, xr);
			row(x,i) = xr;
			//row(x,i) -= prod(LU.get_block(i,2), row(x,i+1));
		}
		b = LU.get_block(i,1);
		column(xrm,0) = row(x,i);
		atlas::lu_substitute_transpose(b, get_pivots(i), xrm);
		row(x,i) = column(xrm,0);
	}

	//Forward substitution pass
	// x[i] = d[i] - anew[i+1] x[i+1]
	for(int i=nrows-2; i>-1; i--)
	{
		xr = row(x,i);
		atlas::gemv(CblasTrans, -1, LU.get_block(i+1,0), row(x,i+1), 1, xr);
		row(x,i) = xr;
		//row(x,i) -= prod(LU.get_block(i+1,0), row(x,i));
	}

};



/*
	Python interface routines
*/
template<class T>
void tridiagonal_block_lu<T>::py_lu(numpy::array& Ain, unsigned nblock, T shift=0)
{
	//Check type
	if( check_type<T>(Ain) ) return;

	//Create wrapped matrix
	T* data_ptr = static_cast<T*>(PyArray_DATA(Ain.ptr()));
	std::vector<unsigned> Asize = shape(Ain);
	
	unsigned nrows = Asize[0]/nblock;
	unsigned ncols = Asize[1]/nblock;
	
	//attach data to block matrix
	LU.assign(nrows, ncols, nblock, data_ptr);

	//Perform LU decomposition
	lu(LU, shift);
};

template<class T>
void tridiagonal_block_lu<T>::py_lu_update(numpy::array& Ain, int uprows)
{
	//Check type
	if( check_type<T>(Ain) ) return;

	//Create wrapped matrix
	T* data_ptr = static_cast<T*>(PyArray_DATA(Ain.ptr()));
	std::vector<unsigned> Asize = shape(Ain);
	
	unsigned nblock = LU.get_nblock();
	unsigned nrows = Asize[0]/nblock;
	unsigned ncols = Asize[1]/nblock;

	//attach data to block matrix
	block_array_wrapped<T> Aupdate(nrows, ncols, nblock, data_ptr);

	//Perform LU decomposition
	lu_update(Aupdate, uprows);
};


template<class T>
void tridiagonal_block_lu<T>::py_solve(numpy::array& x)
{
	typedef typename ublas_extern<T>::adaptor_type ad_t;
	typedef typename ublas_extern<T>::external_matrix_type exmat_t;

	//Check type
	if( check_type<T>(x) ) return;

	//Check size
	unsigned nrows = LU.get_nrows();
	unsigned nblock = LU.get_nblock();
	if( length(x)!=nrows*nblock )
	{
		std::cout << "UBLOCKLU Error: Shape incorrect, A is " << LU.get_nrows() << "x" \
			<< LU.get_nblock() << ", x has " << length(x) << std::endl;
		return;
	}

	T* data_ptr = static_cast<T*>(PyArray_DATA(x.ptr()));
	ad_t shared_adaptor(length(x), data_ptr);
	exmat_t xb(nrows, nblock, shared_adaptor);
	//exmat_t xb = numpy_to_ublas<T>(x);

	solve(xb);
};

template<class T>
void tridiagonal_block_lu<T>::py_solve_transpose(numpy::array& x)
{
	typedef typename ublas_extern<T>::adaptor_type ad_t;
	typedef typename ublas_extern<T>::external_matrix_type exmat_t;

	//Check type
	if( check_type<T>(x) ) return;

	//Check size
	unsigned nrows = LU.get_nrows();
	unsigned nblock = LU.get_nblock();
	if( length(x)!=nrows*nblock )
	{
		std::cout << "UBLOCKLU Error: Shape incorrect, A is " << LU.get_nrows() << "x" \
			<< LU.get_nblock() << ", x has " << length(x) << std::endl;
		return;
	}

	T* data_ptr = static_cast<T*>(PyArray_DATA(x.ptr()));
	ad_t shared_adaptor(length(x), data_ptr);
	exmat_t xb(nrows, nblock, shared_adaptor);
	//exmat_t xb = numpy_to_ublas<T>(x);

	solve_transpose(xb);
};



#endif
