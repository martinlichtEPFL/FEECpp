#ifndef INCLUDEGUARD_DENSE_FACTORIZATION_HPP
#define INCLUDEGUARD_DENSE_FACTORIZATION_HPP

#include "../basic.hpp"

#include "densematrix.hpp"


DenseMatrix GaussJordan( DenseMatrix mat );

DenseMatrix GaussJordanInplace( DenseMatrix mat, bool pivoting = true );


DenseMatrix CholeskyDecomposition( const DenseMatrix& src );

DenseMatrix CholeskyDecompositionBanachchiewicz( const DenseMatrix& src );


void QRFactorization( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R );

void LQFactorization( const DenseMatrix& A, DenseMatrix& L, DenseMatrix& Q );

FloatVector QRIteration( DenseMatrix A, int repetitions = 100 );

FloatVector SolveOverconstrained( const DenseMatrix& A, const FloatVector& v );




 

#endif
