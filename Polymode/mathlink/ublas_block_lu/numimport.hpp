/*
	Conversion routines between numpy and boost arrays
	Currently uses connection to boost-python
*/

#ifndef NUMIMPORT_H
#define NUMIMPORT_H

//Python
#include <Python.h>

//Boost ublas
#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <boost/numeric/ublas/storage.hpp>
#include <boost/numeric/ublas/io.hpp>

//Boost:
#include <boost/python.hpp>
#include <boost/format.hpp>
#include <boost/numeric/ublas/matrix.hpp>

//Numpy:
#include <numpy/noprefix.h>
#include <numpy/arrayobject.h>

//std:
#include <iostream>
#include <vector>

namespace ublas = boost::numeric::ublas;
namespace numpy = boost::python::numeric;

using namespace boost::python;

/* Helper class for soem ublas type names */
template<class T>
struct ublas_extern
{
public:
	typedef typename ublas::shallow_array_adaptor<T> adaptor_type;
	typedef typename ublas::matrix<T, ublas::row_major, adaptor_type> external_matrix_type;
};

//Helper functions to make quick vectors
template<class T>
std::vector<T> zip_vector(T a, T b)
 {
    std::vector<T> vec(2);
	vec[0]=a;
	vec[1]=b;
	return vec;
  };

template<class T>
std::vector<T> zip_vector(T a, T b, T c)
 {
    std::vector<T> vec(3);
	vec[0]=a;
	vec[1]=b;
	vec[2]=b;
	return vec;
  };

//Get an element of the array
template<class T> T get_element(numpy::array& A,int i)
{ return *static_cast<T*>(PyArray_GETPTR1(A.ptr(),i)); }
template<class T> T get_element(numpy::array& A,int i,int j)
{ return *static_cast<T*>(PyArray_GETPTR2(A.ptr(),i,j)); }
template<class T> T get_element(numpy::array& A,int i,int j,int k)
{ return *static_cast<T*>(PyArray_GETPTR3(A.ptr(),i,j,k)); }

//Set an element of the array
template<class T> void set_element(numpy::array& A,int i,T x)
{ *static_cast<T*>(PyArray_GETPTR1(A.ptr(),i)) = x; }
template<class T> void set_element(numpy::array& A,int i,int j,T x)
{ *static_cast<T*>(PyArray_GETPTR2(A.ptr(),i,j)) = x; }
template<class T> void set_element(numpy::array& A,int i,int j,int k,T x)
{ *static_cast<T*>(PyArray_GETPTR3(A.ptr(),i,j,k)) = x; }

//Class to convert template arguments to NUMPY typecodes
//Note this does NOT produce a compile-time error if there is no
// matching definition.
template<class T> struct npy_types
{
	static NPY_TYPES typecode;
	static std::string typedesc;
};

template<class T>
bool check_type(object& obj)
{
  if( PyArray_TYPE(obj.ptr()) != npy_types<T>::typecode )
	 { std::cout << "Array not correct data type, expected a "
		    << npy_types<T>::typedesc << std::endl;
	   return 1;
	 }
  return 0;
};

template<class T>
bool istype(object& obj)
{
  return PyArray_TYPE(obj) == npy_types<T>::typecode;
};


template<class T>
std::vector<T> tuple_to_vec(tuple t)
{
   std::vector<T> vec;
    for (int i = 0; i < len(t); i++)
		vec.push_back( extract<T>(t[i]) );
   return vec;
};

unsigned ndim(numpy::array&);
std::vector<unsigned> shape(numpy::array&);
unsigned length(numeric::array&);

template<class T>
ublas::matrix<T> numpy_to_ublas(numpy::array& Ain)
{
	typedef typename ublas_extern<T>::adaptor_type ad_t;
	typedef typename ublas_extern<T>::external_matrix_type exmat_t;

	check_type<T>(Ain);

	T* data_ptr = static_cast<T*>(PyArray_DATA(Ain.ptr()));
	std::vector<unsigned> Asize = shape(Ain);
	
	ad_t shared_adaptor(Asize[0]*Asize[1], data_ptr);
	exmat_t Aout(Asize[0], Asize[1], shared_adaptor);
	
	return Aout;
};

template<class T>
numpy::array ublas_to_numpy(const ublas::matrix<T>& Ain)
{
	T* data_ptr = const_cast<T*>(Ain.data());
	double data_shape[2];
	data_shape[0] = Ain.size1();
	data_shape[1] = Ain.size2();

	object obj(handle<>(PyArray_SimpleNewFromData(2, &data_shape, npy_types<T>::typecode, data_ptr)));
   
   return extract<numpy::array>(obj);
};


#endif
