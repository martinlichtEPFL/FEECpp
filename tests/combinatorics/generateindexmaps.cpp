

/**/

#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"
#include "../../combinatorics/indexmap.hpp"
#include "../../combinatorics/generateindexmaps.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Index Map Generators" << nl;
    
    if(true)
    {
            
        LOG << "1. Testing generator for empty index map" << nl;
        
        const IndexRange from( 1,-3 );
        const IndexRange targ( 2,5 );
        
        const std::vector<IndexMap>  all = generateEmptyMap( from, targ );
        
        for( const IndexMap& im : all ) LOG << im << nl;

        assert( all.size() == 1 );
        
        LOG << "Tested" << nl;
        
    }
    
    if(true)
    {
        
        LOG << "2. Testing generator for all: [1,3] -> [1,3]" << nl;
        
        IndexRange bereich( 1,3 );
        std::vector<IndexMap> all = generateIndexMaps( bereich, bereich );

        for( const IndexMap& im : all ) LOG << im << nl;

        assert( all.size() == 27 );
        
        LOG << "Tested" << nl;
        
    }
    
    if(true)
    {
        
        LOG << "3. Testing generator for permutations of [1,3]" << nl;
        
        IndexRange bereich( 1,3 );
        std::vector<IndexMap>  all = generatePermutations( bereich );

        for( const IndexMap& im : all ) LOG << im << nl << "signum: " << signPermutation( im ) << nl;
        
        assert( all.size() == 6 );
        
        LOG << "Tested" << nl;
        
    }
    
    if(true)
    {
        
        LOG << "4. Testing generator for sigmas [1,3] -> [2,5]" << nl;
        
        IndexRange from( 1,3 );
        IndexRange targ( 2,5 );
        std::vector<IndexMap>  all = generateSigmas( from, targ );

        for( const IndexMap& im : all ) LOG << im << nl;
        
        assert( all.size() == 4 );
        
        LOG << "Tested" << nl;
            
    }
    
    
    
    
    if(true)
    {
            
        LOG << "5. Testing generator for general index mappings" << nl;
        
        const std::vector<int> Ns = { -1, 0, 1, 2, 3, 4 };
        const std::vector<int> Ks = { -1, 0, 1, 2, 3, 4 };
        
        
        for( const int& N : Ns ) 
        for( const int& K : Ks ) 
        {
        
            const IndexRange source( 1, N );
            const IndexRange target( 1, K );
            
            const std::vector<IndexMap> all = generateIndexMaps( source, target );
            
            if( source == target )
                assert( all == generateIndexMaps( source) );

            if( source.cardinality() == 0 and target.cardinality() == 0 )
                assert( all.size() == 1 );
            else
                assert( all.size() == power_integer( target.cardinality(), source.cardinality() ) );
            
            for( int i = 0; i < all.size(); i++ )
            for( int j = 0; j < all.size(); j++ )
                if( i != j )
                    assert( all[i] != all[j] );
            
            for( const IndexMap& im : all ) LOG << im << nl;
            
        }
        
        LOG << "Tested" << nl;
        
    }
    
    if(true)
    {
            
        LOG << "6. Testing generator for permutations" << nl;
        
        const std::vector<int> Ns = { -1, 0, 1, 2, 3, 5 };
        
        
        for( const int& N : Ns ) {
            
            const IndexRange bereich( 1, N );
            const std::vector<IndexMap>  all = generatePermutations( bereich );

            // check cardinality 
            if( bereich.cardinality() > 0 )
                assert( all.size() > 0 );
            
            assert( all.size() == factorial_integer( bereich.cardinality() ) );
            
            // check they are all different 
            for( int i = 0; i < all.size(); i++ )
            for( int j = 0; j < all.size(); j++ )
                if( i != j )
                    assert( all[i] != all[j] );
            
            // print their signa
            for( const IndexMap& im : all ) LOG << im << nl << "signum: " << signPermutation( im ) << nl;
            
            // count the even permutations 
            int number_of_even_permutations = 0;
            
            for( const IndexMap& im : all ) {
                
                assert( im.getSourceRange() == im.getTargetRange() );
                
                for( const int i : bereich ) {
                
                    assert( bereich.min() <= im[i] and im[i] <= bereich.max() );
                
                    for( const int j : bereich  )
                        if( i != j )
                            assert( im[i] != im[j] );
                
                }
                
                assert( im.isbijective() );
                
                if( signPermutation( im ) == 1 ) 
                    number_of_even_permutations++;
            }
            
            if( bereich.cardinality() > 1 )
                assert( number_of_even_permutations * 2 == all.size() );
            else
                assert( number_of_even_permutations     == all.size() );
            
        }
        
        LOG << "Tested" << nl;
        
    }
    
    if(true)
    {
            
        LOG << "Testing generator for sigmas [1,N] -> [2,L]" << nl;
        
        const std::vector<int> Ns = { -1, 0, 1, 2, 3, 5 };
        const std::vector<int> Ls = {  1, 2, 3, 4, 5, 7, 9 };
        
        
        for( const int& N : Ns ) 
        for( const int& L : Ls ) 
        {
            
            IndexRange source( 1, N );
            IndexRange target( 2, L );
            
            std::vector<IndexMap>  all = generateSigmas( source, target );
            
            if( source.cardinality() > target.cardinality() )
                assert( all.size() == 0 );
            
            assert( all.size() == binomial_integer( target.cardinality(), source.cardinality() ) );
            
            for( const IndexMap& sigma : all ) {
            
                assert( sigma.getSourceRange() == source );
                assert( sigma.getTargetRange() == target );
            
                for( const int i : source ) {
            
                    assert( target.min() <= sigma[i] and sigma[i] <= target.max() );
            
                    if( i != source.max() )
                        assert( sigma[i] < sigma[i+1] );
            
                }
            
                assert( sigma.isstrictlyascending() );
                
            }
            
            for( const IndexMap& sigma : all ) LOG << sigma << nl;
            
        }
        
        LOG << "Tested" << nl;
            
    }
    
    if(true)
    {
        
        LOG << "Testing generator for sigmas in 2D:" << nl;
        
        const std::vector<IndexRange> sources = { IndexRange(1,1), IndexRange(1,2) };
        const IndexRange target( 0,2 );
        
        for( auto source : sources ) 
        {
            std::vector<IndexMap> sigmas = generateSigmas( source, target );
            for( const IndexMap& im : sigmas ) LOG << im << nl;
            assert( sigmas.size() == binomial_integer( 3, source.max() ) );
        }
        LOG << "Tested" << nl;
        
    }
    
    if(true)
    {
        
        LOG << "Testing generator for sigmas in 3D:" << nl;
        
        const std::vector<IndexRange> sources = { IndexRange(1,1), IndexRange(1,2), IndexRange(1,3) };
        const IndexRange target( 0,3 );
        
        for( auto source : sources ) 
        {
            std::vector<IndexMap> sigmas = generateSigmas( source, target );
            for( const IndexMap& im : sigmas ) LOG << im << nl;
            assert( sigmas.size() == binomial_integer( 4, source.max() ) );
        }
        LOG << "Tested" << nl;
        
    }
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
