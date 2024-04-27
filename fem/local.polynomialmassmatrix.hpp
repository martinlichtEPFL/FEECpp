#ifndef INCLUDEGUARD_FEM_POLYNOMIALMASSMATRIX_HPP
#define INCLUDEGUARD_FEM_POLYNOMIALMASSMATRIX_HPP

#include "../combinatorics/multiindex.hpp"
#include "../dense/densematrix.hpp"



DenseMatrix polynomialmassmatrix( int n, int r );

DenseMatrix polynomialmassmatrix( int n, int r, const MultiIndex& base );

std::vector<DenseMatrix> polynomialmassmatrices_per_multiindex( int n, int r, int s );

std::vector<DenseMatrix> polynomialmassmatrices_per_lagrangepoint( int n, int r, int s );


#endif
