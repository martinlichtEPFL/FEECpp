#include <algorithm>
#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"


#include "../fem/indexfunctions.hpp"

std::vector< std::pair<MultiIndex,IndexMap> > ListOfSullivanIndices( int n, int k, int r )
{
    
    // check whether the parameters are right 
    
    assert( n >= 0 );
    assert( k >= 0 );
    assert( r >= 0 );
    assert( r >= 1 or k >= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    // if k > n, then there is nothing to show 
    if( k > n ) 
        return std::vector< std::pair<MultiIndex,IndexMap> >();
    
    // if the polynomial degree is zero, check we do volume forms 
    if( r == 0 ) {
        assert( k >= n );
        MultiIndex foo = ZeroMultiIndex( IndexRange(0,n) );
        IndexMap   bar = IndexMap( IndexRange(1,n), IndexRange(0,n), [n](int i) -> int { assert( 1 <= i and i <= n ); return i; } );
        std::vector< std::pair<MultiIndex,IndexMap> > ret { std::pair<MultiIndex,IndexMap>( foo, bar ) };
        return ret;
    }
    
    
    // Auxiliary calculations and preparations
    
    // List of the numbers 0..n
    const std::vector<int> N = [&n]()->std::vector<int>{ 
        std::vector<int> ret(n+1); 
        for( int i = 0; i <= n; i++ ) ret[i] = i;
        return ret;
    }();
    
    const std::vector<MultiIndex> alphas = generateMultiIndices( IndexRange( 0, n ), r );
    
    const std::vector<IndexMap>   sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n ) );
    
    std::vector< std::pair<MultiIndex,IndexMap> > ret;
    
    //  [ size of set taken from Acta paper : Theorem 4.10 ]
    int computed_length = binomial_integer( r-1, n-k ) * binomial_integer( r+k, k );
    
    
    // filter out the basis forms 
    for( const MultiIndex& alpha : alphas )
    for( const IndexMap&   sigma : sigmas )
    {
        
        // First, check that every p in 0..n is contained in the ranges of alpha and/or sigma
        bool b1 = std::all_of( N.begin(), N.end(), 
                               [ &alpha, &sigma ]( int p ) -> bool { 
                                   return sigma.has_value_in_range(p) or alpha[p] > 0;
                                   }
                             );
        
        // Second, check that min[alpha] not in [sigma]
        bool b2 = not sigma.has_value_in_range( alpha.min() );
        
        // if both criteria are satisfied, then save that one
        if( b1 and b2 )
            ret.push_back( std::pair<MultiIndex,IndexMap>( alpha, sigma ) );
        
    }
    
    assert( ret.size() == computed_length );
    
    return ret;
    
}



std::vector< std::pair<MultiIndex,IndexMap> > ListOfWhitneyIndices( int n, int k, int r )
{
    
    // check whether the parameters are right 
    
    assert( n >= 0 );
    assert( k >= 0 );
    assert( r >= 1 );
    
    // if k > n, then there is nothing to show 
    if( k > n ) 
        return std::vector< std::pair<MultiIndex,IndexMap> >();
    
    // Auxiliary calculations and preparations
    
    // List of the numbers 0..n
    const std::vector<int> N = [&n]()->std::vector<int>{ 
        std::vector<int> ret(n+1); 
        for( int i = 0; i <= n; i++ ) ret[i] = i;
        return ret;
    }();
    
    const std::vector<MultiIndex> alphas = generateMultiIndices( IndexRange( 0, n ), r-1 );
    
    const std::vector<IndexMap>   rhos   = generateSigmas( IndexRange( 0, k ), IndexRange( 0, n ) );
    
    std::vector< std::pair<MultiIndex,IndexMap> > ret;
    
    //  [ size of set taken from Acta paper: Theorem 4.14 ]
    int computed_length = binomial_integer( n, k ) * binomial_integer( r+k-1, n );
    
    
    // filter out the basis forms 
    for( const MultiIndex& alpha : alphas )
    for( const IndexMap&     rho :   rhos )
    {
        
        // First, check that every p in 0..n is contained in the ranges of alpha and/or sigma
        bool b1 = std::all_of( N.begin(), N.end(), 
                               [ &alpha, &rho ]( int p ) -> bool { 
                                   return rho.has_value_in_range(p) or alpha[p] > 0;
                                   }
                             );
        
        // Second, check that min[sigma] = 0
        bool b2 = rho.has_value_in_range( 0 );
        
        // if both criteria are satisfied, then save that one
        if( b1 and b2 )
            ret.push_back( std::pair<MultiIndex,IndexMap>( alpha, rho ) );
        
    }
    
    assert( ret.size() == computed_length );
    
    return ret;
    
}
