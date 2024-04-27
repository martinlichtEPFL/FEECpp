
#include "densematrix.hpp"
#include "functions.hpp"
#include "simplesolver.hpp"
#include "factorization.hpp"

#include <new>

 // LR factorization, column pivot, 
 
DenseMatrix GaussJordan( DenseMatrix mat )
{
    
    assert( mat.issquare() );
    
    const int n = mat.getdimout();
//     std::vector<int> pivotcol( n, -17 );
//     std::vector<int> colperm ( n, -17 );
//     for( int c = 0; c < n; c++ ) colperm[c] = c;

    
    DenseMatrix ret(n);
    ret.unitmatrix();
    
    // 1. eliminate lower triangular part, save coeffecients
    
    for( int i = 0; i < n; i++ ) {
        
        for( int k = 0; k < n; k++ ) { // each
            
            if( i == k ) continue; 
            
            Float coeff = - mat( k, i ) / mat( i, i );
            
            for( int j = i; j < n; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
            }

            for( int j = 0; j <= i; j++ ) {
                ret( k, j ) = ret( k, j ) + coeff * ret( i, j );
            }
            
        }
        
        Float coeff = 1. / mat(i,i);
        
        for( int j = 0; j <= i; j++ ) {
            ret(i,j) *= coeff;
        }
        
        for( int j = i; j < n; j++ ) {
            mat(i,j) *= coeff;
        }
            
    }
    
// //     2. normalize the diagonals
//     
//     for( int i = 0; i < n; i++ ) {
//         
//         Float coeff = 1. / mat(i,i);
//         
//         for( int k = 0; k < n; k++ ) {
//             mat(i,k) *= coeff;
//             ret(i,k) *= coeff;
//         }
        
//         ret(i,i) = coeff;
//     }
    
    // finished!
    
    // LOG << mat;
    // LOG << ret;
    
    return ret;
}

 

 
 
 
 
DenseMatrix GaussJordanInplace( DenseMatrix mat, bool pivoting )
{
    
    assert( mat.issquare() );
    
    const int n = mat.getdimout();
    
    int* pivots = nullptr;
    if(pivoting) pivots = new (std::nothrow) int[n];
    
    for( int i = 0; i < n; i++ ) {
        
        if( pivoting ) {
            
            int c_max = i;
            for( int c = i+1; c < n; c++ )
                if( absolute(mat(i,c)) > absolute(mat(i,c_max)) ) 
                    c_max = c;
            
            pivots[i] = c_max;
            mat.swapcolumn( c_max, i );
            
        }
        
        for( int k = 0; k < n; k++ ) { // each
            
            if( i == k ) continue; 
            
            assert( absolute(mat(i,i)) != 0.0 );
            
            Float coeff = - mat( k, i ) / mat( i, i );
            
            for( int j = i+1; j < n; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
            }

            mat( k, i ) = coeff;
            
            for( int j = 0; j < i; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
            }
            
        }
        
        Float coeff = 1. / mat(i,i);
        
        for( int j = 0; j < i; j++ ) {
            mat(i,j) *= coeff;
        }
        
        mat(i,i) = coeff;
        
        for( int j = i+1; j < n; j++ ) {
            mat(i,j) *= coeff;
        }
            
    }
    
    if( pivoting ) {
        for( int i = n-1; i >= 0; i-- )
//         for( int i = 0; i < n; i++ ) 
        {
//             LOG << "swap " << i << space << pivots[i] << nl;
            mat.swaprow( i, pivots[i] );
        }
    }
    
    // finished!
    
    if( pivoting ) delete[] pivots;
    
    return mat;
}








DenseMatrix CholeskyDecomposition( const DenseMatrix& A )
{
    return CholeskyDecompositionBanachchiewicz( A );
}


DenseMatrix CholeskyDecompositionBanachchiewicz( const DenseMatrix& A )
{
    A.check();
    assert( A.issquare() );

    DenseMatrix L = A;
    L.set( 0. );
    const int dim = A.getdimout();

    for( int r = 0; r < dim; r++ ){
        
        for( int c = 0; c < r; c++ ) {
        
            L(r,c) = A(r,c);
            
            for( int k = 0; k < c; k++ )
                L(r,c) -= L(r,k) * L(c,k);
            
            L(r,c) /= L(c,c);
            
            assert( absolute( L(c,c) ) > 0. );
        }
        
        L(r,r) = A(r,r);
        
        for( int k = 0; k < r; k++ )
            L(r,r) -= L(r,k) * L(r,k);
        
        assert( L(r,r) >= 0. );
        
        L(r,r) = std::sqrt( L(r,r) );
        
    }

    L.check();
    return L;
}




void QRFactorization( const DenseMatrix& A, DenseMatrix& Q, DenseMatrix& R )
{
    A.check();
    Q.check();
    R.check();
    assert( A.getdimout() == Q.getdimout() );
    assert( Q.getdimin()  == R.getdimout() );
    assert( A.getdimin()  == R.getdimin()  ); 
    assert( A.getdimin()  <= A.getdimout() );
    assert( R.issquare() );
    
    R.zeromatrix();
    
    for( int c = 0; c < A.getdimin(); c++ ) {
        FloatVector u = A.getcolumn(c);
        for( int j = 0; j < c; j++ ){
                R(j,c) = u * Q.getcolumn(j);
                u -= R(j,c) * Q.getcolumn(j);
        }
        R(c,c) = std::sqrt( u*u );
        Q.setcolumn( c, u / R(c,c) );
    }
    
}

void LQFactorization( const DenseMatrix& A, DenseMatrix& L, DenseMatrix& Q )
{
    A.check();
    L.check();
    Q.check();
    
    assert( A.getdimout() == L.getdimout() );
    assert( A.getdimin()  == Q.getdimin()  ); 
    assert( Q.getdimout() == L.getdimin() );
    assert( A.getdimout() <= A.getdimin() );
    assert( L.issquare() );
    
    L.zeromatrix();
    
    for( int r = 0; r < A.getdimout(); r++ ) {
        FloatVector v = A.getrow(r);
        for( int i = 0; i < r; i++ ){
                L(r,i) = v * Q.getrow(i);
                v -= L(r,i) * Q.getrow(i);
        }
        L(r,r) = std::sqrt( v*v );
        Q.setrow( r, v / L(r,r) );
    }
    
}


FloatVector SolveOverconstrained( const DenseMatrix& A, const FloatVector& b )
{
    assert( A.getdimout() == b.getdimension() );

    DenseMatrix Q( A.getdimout(), A.getdimin() );
    DenseMatrix R( A.getdimin() );
    QRFactorization(A,Q,R);

    FloatVector x( A.getdimin() );

    UpperTriangularSolve( R, x, Transpose(Q) * b );

    return x; 
}


FloatVector QRIteration( DenseMatrix A, int repetitions ) 
{
    assert( A.issquare() && A.issymmetric() );
    
    const int dim = A.getdimin();
    
    while( repetitions --> 0 )
    {
        DenseMatrix Q(dim,dim), R(dim,dim);        
        QRFactorization( A, Q, R ); A = R * Q;
    }

    return A.getDiagonal();
}




