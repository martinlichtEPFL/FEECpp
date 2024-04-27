
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../dense/densematrix.hpp"
#include "../fem/utilities.hpp"

#include "../fem/local.polynomialmassmatrix.hpp"



DenseMatrix polynomialmassmatrix( int n, int r )
{
    assert( n >= 0 && r >= 0 );
    
    // create the multiindices and the return value matrix
    
    std::vector<MultiIndex> multis = generateMultiIndices( IndexRange(0,n), r );
    
    const int N = multis.size();
    
    assert( N == binomial_integer( n + r , r ) );
    
    DenseMatrix ret( N, N, notanumber );
    
    // create the entries  
    
    for( int i = 0; i < N; i++ )
    for( int j = 0; j < N; j++ )
    {
        
        MultiIndex alpha = multis[i] + multis[j];
        
        //  alpha! / (n+|alpha|)!
        assert( absolute( alpha ) == 2*r );
        
        ret( i, j ) = factorial_numerical(n) * alpha.factorial_numerical() / factorial_numerical( n + 2*r ); 
        
        if( ret( i, j ) <= 0. ) {
            LOG << multis[i] << multis[j] << ret (i, j );
            unreachable();
        }

    }
    
    assert( ret.isfinite() );
    
    if( not ret.ispositive() )
        LOG << ret;
    
    assert( ret.isnonnegative() );
    assert( ret.ispositive() );
    
    return ret;
}






DenseMatrix polynomialmassmatrix( int n, int r, const MultiIndex& base )
{
    assert( n >= 0 && r >= 0 );
    
    // create the multiindices and the return value matrix
    
    std::vector<MultiIndex> multis = generateMultiIndices( IndexRange(0,n), r );
    
    const int N = multis.size();
    
    assert( N == binomial_integer( n + r , r ) );
    
    DenseMatrix ret( N, N, notanumber );
    
    // create the entries  
    
    for( int i = 0; i < N; i++ )
    for( int j = 0; j < N; j++ )
    {
        
        MultiIndex alpha = multis[i] + multis[j] + base;
        
        //  alpha! / (n+|alpha|)!
        assert( absolute( alpha ) == 2*r + absolute(base) );
        
        ret( i, j ) = factorial_numerical(n) * alpha.factorial_numerical()
                      / factorial_numerical( n + 2*r + absolute(base) ); 
        
        if( ret( i, j ) <= 0. ) {
            LOG << multis[i] << multis[j] << ret (i, j );
            unreachable();
        }

    }
    
    assert( ret.isfinite() );
    
    if( not ret.ispositive() )
        LOG << ret;
    
    assert( ret.isnonnegative() );
    assert( ret.ispositive() );
    
    return ret;
}





std::vector<DenseMatrix> polynomialmassmatrices_per_multiindex( int n, int r, int s )
{
    assert( n >= 0 && r >= 0 && s >= 0 );
    
    // create the multiindices and the return value matrix
    
    std::vector<MultiIndex> multis = generateMultiIndices( IndexRange(0,n), s );

    const int N = multis.size();

    assert( N == binomial_integer( n + s , s ) );
    
    std::vector<DenseMatrix> ret( N, DenseMatrix( binomial_integer( n + r , r ) ) );

    for( int i = 0; i < N; i++ )
        ret[i] = polynomialmassmatrix( n, r, multis[i] );
    
    return ret;
}


std::vector<DenseMatrix> polynomialmassmatrices_per_lagrangepoint( int n, int r, int s )
{
    assert( n >= 0 && r >= 0 && s >= 0 );

    int N = binomial_integer( n+s, n );
    
    std::vector<DenseMatrix> massmatrices = polynomialmassmatrices_per_multiindex( n, r, s );
    
    DenseMatrix coefficients = LagrangePolynomialCoefficients( n, s );
    // each column contains the coefficients for the monomials

    std::vector<DenseMatrix> ret( N, DenseMatrix( binomial_integer( n + r , r ), (Float)0. ) );

    for( int i = 0; i < N; i++ ) // i -> index of lagrangepoint 
    for( int j = 0; j < N; j++ ) // j -> index of monomial in lagrange polynomial
        ret[i] += coefficients(j,i) * massmatrices[j];
    
    return ret;
}

