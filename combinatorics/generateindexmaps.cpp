
#include "generateindexmaps.hpp"

#include <vector>

#include <algorithm>
#include <random>

#include "../basic.hpp"
#include "indexrange.hpp"
#include "indexmap.hpp"




std::vector<IndexMap> generateEmptyMap( const IndexRange& from, const IndexRange& to )
{
    
    from.check();
    to.check();
    assert( from.isempty() );
    
    const IndexMap myempty = IndexMap( from, to, std::vector<int>(0) );
    myempty.check();
    
    std::vector<IndexMap> ret;
    ret.push_back( myempty );

    if( /* DISABLES CODE */ ( false ) ){ 
        std::mt19937 g( 234567891 );
        std::shuffle( ret.begin(), ret.end(), g );
    }
    
    return ret;    
}





std::vector<IndexMap> generateIndexMaps( const IndexRange& range )
{
    return generateIndexMaps( range, range );
}

std::vector<IndexMap> generateIndexMaps( const IndexRange& from, const IndexRange& to )
{
    
    from.check();
    to.check();
    
    if( from.isempty() )
        return generateEmptyMap( from, to );
    
    if( to.isempty() )
        return std::vector<IndexMap>();
    
    assert( from.cardinality() > 0 );
    assert( to.cardinality() > 0 );
    
    const int num = power_integer( to.cardinality(), from.cardinality() );
    
    std::vector<IndexMap> ret;
    ret.reserve( num );
    
    for( int it = 0; it < num; it++ )
    {
        
        std::vector<int> values( from.cardinality() );
        int it_clone = it;
        
        for( int digit = 0; digit < from.cardinality(); digit++ )
        {
            values[digit] = to.min() + it_clone % to.cardinality();
            it_clone /= to.cardinality();
        }

        ret.emplace_back( from, to, values );
        
    }
    
    for( const auto& foo : ret )
        foo.check();
    
    if( /* DISABLES CODE */ ( false ) ){ 
        std::mt19937 g( 345678912 );
        std::shuffle( ret.begin(), ret.end(), g );
    }
    
    return ret;
    
}











std::vector<IndexMap> generatePermutations( const IndexRange& ir )
{
    
    ir.check();
    const std::vector<IndexMap> allmappings = generateIndexMaps( ir, ir );
    
    std::vector<IndexMap> ret;
    ret.reserve( factorial_integer( ir.cardinality() ) );
    
    for( const auto& mapping : allmappings )
        if( mapping.isbijective() )
        ret.push_back( mapping );
            
    for( const auto& perm : ret )
        assert( (perm.check(),true) && perm.isbijective() );
    
    if( /* DISABLES CODE */ ( false ) ){ 
        std::mt19937 g( 456789123 );
        std::shuffle( ret.begin(), ret.end(), g );
    }
    
    return ret;
    
}

int signPermutation( const IndexMap& im )
{
    
    im.check();
    assert( im.isbijective() );
    assert( im.getSourceRange() == im.getTargetRange() );
    
    const IndexRange& ir = im.getSourceRange();
    
    bool is_even = true;
    
//     for( int s = ir.min(); s <= ir.max(); s++ )
//     for( int t = s+1;      t <= ir.max(); t++ )
//     {
//         assert( s < t );
//         assert( im[s] != im[t] );
//         if( im[s] > im[t] ) is_even = not is_even; 
//     }

    for( const int s : ir )
    for( const int t : ir )
    {
        if( s >= t ) continue;
        assert( s < t );
        assert( im[s] != im[t] );
        if( im[s] > im[t] )
            is_even = not is_even; 
    }
    
    int ret = is_even ? 1 : -1;
    
    
//     int zaehler = 1;
//     int nenner = 1;
//     
//     for( int s = ir.min(); s <= ir.max(); s++ )
//     for( int t = s+1; t <= ir.max(); t++ )
//     {
//         nenner *= ( t - s );
//         zaehler *= ( im[ t ] - im[ s ] );
//     }
// 
//     assert( zaehler % nenner == 0 );
//     int new_ret = zaehler / nenner;
//     assert( new_ret == 1 || new_ret == -1 );
//     
//     assert( new_ret == ret );
    
    
    return ret;
    
}










std::vector<IndexMap> generateSigmas( const IndexRange& from, const IndexRange& to )
{
    
    from.check();
    to.check();
    
    const std::vector<IndexMap> allmappings = generateIndexMaps( from, to );
    
    std::vector<IndexMap> ret;
    
    if( from.isempty() ) 
        assert( allmappings.size() == 1 ); 
    
    ret.reserve( binomial_integer( to.cardinality(), from.cardinality() ) );
    
    for( const auto& mapping : allmappings )
        if( mapping.isstrictlyascending() )
            ret.push_back( mapping );
            
    for( const auto& sigma : ret ) {
        sigma.check();
        assert( sigma.isstrictlyascending() );
    }
    
    if( /* DISABLES CODE */ ( false ) ){ 
        std::mt19937 g( 567891234 );
        std::shuffle( ret.begin(), ret.end(), g );
    }
    
    return ret;
    
}



