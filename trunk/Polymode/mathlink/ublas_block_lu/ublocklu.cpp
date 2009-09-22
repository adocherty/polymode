/*
    Boost python interface to UBlas Block LU c++ module
*/

#include "ublocklu.hpp"
#include "numimport.hpp"

using namespace boost::python;

BOOST_PYTHON_MODULE(ublocklu)
{
    class_<tridiagonal_block_lu<double> >("dblocklu", init< optional<bool> >())
       .def("__call__",&tridiagonal_block_lu<double>::py_lu)
       .def("update",&tridiagonal_block_lu<double>::py_lu_update)
       .def("solve",&tridiagonal_block_lu<double>::py_solve)
       .def("solve_transpose",&tridiagonal_block_lu<double>::py_solve_transpose)
        ;

    class_<tridiagonal_block_lu<complex_t> >("cblocklu", init< optional<bool> >())
       .def("__call__",&tridiagonal_block_lu<complex_t>::py_lu)
       .def("update",&tridiagonal_block_lu<complex_t>::py_lu_update)
       .def("solve",&tridiagonal_block_lu<complex_t>::py_solve)
       .def("solve_transpose",&tridiagonal_block_lu<complex_t>::py_solve_transpose)
        ;

    //Important! .. Direct numpy calls won't work if we don't have this.
    import_array();
    numeric::array::set_module_and_type("numpy", "ndarray");

}
