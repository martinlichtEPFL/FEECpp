#ifndef INCLUDEGUARD_DENSE_FUNCTIONS_HPP
#define INCLUDEGUARD_DENSE_FUNCTIONS_HPP

#include "../basic.hpp"

#include "densematrix.hpp"

 
 // Tranpose, giving new matrix 
 // Tranpose, in-place

DenseMatrix Transpose( const DenseMatrix& ); 

// void TransposeInSitu( DenseMatrix& ); // The interface does not allow changing the dimension of the matrix

DenseMatrix TransposeSquare( const DenseMatrix& ); 

void TransposeSquareInSitu( DenseMatrix& ); 




// remove single rows or columns
// DenseMatrix skip_row( int i, const DenseMatrix& mat );
// DenseMatrix skip_column( int i, const DenseMatrix& mat );






 // determinant 
 // cofactor matrix 
 // Inverse 
 // Inverse and determinant 
 // SubdeterminantMatrix

Float Determinant( const DenseMatrix& );

DenseMatrix CofactorMatrix( const DenseMatrix& ); 

DenseMatrix Inverse( DenseMatrix );

void Inverse_InSitu( DenseMatrix& );

// void InverseAndDeterminant( const DenseMatrix&, DenseMatrix&, Float& );

DenseMatrix SubdeterminantMatrixSquare( const DenseMatrix& A, int k );

DenseMatrix SubdeterminantMatrix( const DenseMatrix& A, int k );



Float Determinant_laplaceexpansion( const DenseMatrix& );

Float Determinant_gauss( DenseMatrix );

Float Determinant_bareiss( DenseMatrix A );

void Inverse_CramersRule_InSitu( DenseMatrix& );

void Inverse_gauss_InSitu( DenseMatrix&, bool pivoting = true );


DenseMatrix MatrixTensorProduct( const DenseMatrix& left, const DenseMatrix& right );








#endif
