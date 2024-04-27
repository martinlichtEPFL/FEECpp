#ifndef INCLUDEGUARD_DENSE_SIMPLESOLVER_HPP
#define INCLUDEGUARD_DENSE_SIMPLESOLVER_HPP

#include "../basic.hpp"
#include "../operators/floatvector.hpp"

#include "densematrix.hpp"

 

DenseMatrix DiagonalPart( const DenseMatrix& A );
DenseMatrix DiagonalInverse( const DenseMatrix& A );
void InvertDiagonal( DenseMatrix& A );
Float DiagonalDeterminant( const DenseMatrix& A );
void DiagonalSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );

DenseMatrix LowerTriangularPart( const DenseMatrix& A );
DenseMatrix LowerTriangularInverse( const DenseMatrix& A );
void InvertLowerTriangular( DenseMatrix& A );
Float LowerTriangularDeterminant( const DenseMatrix& A );
void LowerTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );

DenseMatrix UpperTriangularPart( const DenseMatrix& A );
DenseMatrix UpperTriangularInverse( const DenseMatrix& A );
void InvertUpperTriangular( DenseMatrix& A );
Float UpperTriangularDeterminant( const DenseMatrix& A );
void UpperTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );

DenseMatrix LowerUnitTriangularPart( const DenseMatrix& A );
DenseMatrix LowerUnitTriangularInverse( const DenseMatrix& A, bool writediagonalones = false );
void InvertLowerUnitTriangular( DenseMatrix& A, bool writediagonalones = false );
void LowerUnitTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );

DenseMatrix UpperUnitTriangularPart( const DenseMatrix& A );
DenseMatrix UpperUnitTriangularInverse( const DenseMatrix& A, bool writediagonalones = false );
void InvertUpperUnitTriangular( DenseMatrix& A, bool writediagonalones = false );
void UpperUnitTriangularSolve( const DenseMatrix& A, FloatVector& x, const FloatVector& b );



#endif
