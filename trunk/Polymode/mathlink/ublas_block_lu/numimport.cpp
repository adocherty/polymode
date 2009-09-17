//Explicit specializations for the npy_types template:
#include "numimport.hpp"

/*
from numpy/ndarrayobject.h
enum NPY_TYPES {    NPY_BOOL=0,
                    NPY_BYTE, NPY_UBYTE,
                    NPY_SHORT, NPY_USHORT,
                    NPY_INT, NPY_UINT,
                    NPY_LONG, NPY_ULONG,
                    NPY_LONGLONG, NPY_ULONGLONG,
                    NPY_FLOAT, NPY_DOUBLE, NPY_LONGDOUBLE,
                    NPY_CFLOAT, NPY_CDOUBLE, NPY_CLONGDOUBLE,
                    NPY_OBJECT=17,
                    NPY_STRING, NPY_UNICODE,
                    NPY_VOID,
                    NPY_NTYPES,
                    NPY_NOTYPE,
                    NPY_CHAR,
                    NPY_USERDEF=256
};
*/

/* Mappings from NUMPY types to C++ types */
template<> NPY_TYPES npy_types<int>::typecode = NPY_INT;
template<> std::string npy_types<int>::typedesc = "int";

template<> NPY_TYPES npy_types<double>::typecode = NPY_DOUBLE;
template<> std::string npy_types<double>::typedesc = "double";

template<> NPY_TYPES npy_types<float>::typecode = NPY_FLOAT;
template<> std::string npy_types<float>::typedesc = "float";

template<> NPY_TYPES npy_types<std::complex<double> >::typecode = NPY_CDOUBLE;
template<> std::string npy_types<std::complex<double> >::typedesc = "complex double";
template<> NPY_TYPES npy_types<std::complex<float> >::typecode = NPY_CFLOAT;
template<> std::string npy_types<std::complex<float> >::typedesc = "complex float";

//Number of dimentions
unsigned ndim(numeric::array& arr)
{
   return PyArray_NDIM(arr.ptr());
}

//Return shape of array in vector
std::vector<unsigned> shape(numeric::array& arr)
{
  std::vector<unsigned> out_dims;

  intp* dims_ptr = PyArray_DIMS(arr.ptr());

  for (unsigned i = 0; i < ndim(arr); i++){
	unsigned dim = static_cast<unsigned>(*(dims_ptr + i));
    out_dims.push_back(dim);
  }
  return out_dims;
}

//Return total number of elements in array over all dimensions
unsigned length(numeric::array& arr)
{
  intp* dims_ptr = PyArray_DIMS(arr.ptr());
  unsigned length = 1;
  for (unsigned i = 0; i < ndim(arr); i++){
	unsigned dim = static_cast<unsigned>(*(dims_ptr + i));
    length = length*dim;
  }
  return length;
}


