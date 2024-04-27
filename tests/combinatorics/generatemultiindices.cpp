
#include "../../basic.hpp"
#include "../../combinatorics/generatemultiindices.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Multiindex Generators" << nl;
    
    if(true)
    {
            
        LOG << "1. Selected combinations" << nl;
        
        const std::vector<int> mins = { 1, 0, 0, 0, 0, 2,  2 };
        const std::vector<int> maxs = { 4, 2, 2, 2, 2, 2, -3 };
        const std::vector<int> pows = { 2, 2, 3, 0, 1, 7,  3 };
    
        assert( mins.size() == maxs.size() && mins.size() == pows.size() );
    
        for( int m = 0; m < mins.size(); m++ )
        {
            
            int min = mins[m]; int max = maxs[m]; int pow = pows[m];
            
            LOG << "Check: " << min << space << max << space << pow << nl;
        
            IndexRange ir( min, max );
            std::vector<MultiIndex> vmi = generateMultiIndices( ir, pow );

            if( ir.cardinality() == 0 ) {

                LOG << vmi.size() << nl;
                
                for( const auto& mi : vmi ) 
                    LOG << mi << space << mi.absolute() << nl;
                
                Assert( vmi.size() == 0 );
                
                continue;
            }

            assert( vmi.size() == binomial_integer(ir.cardinality() - 1 + pow, ir.cardinality() - 1 ) );
            assert( vmi.size() == binomial_integer(ir.cardinality() - 1 + pow, pow                  ) );
            
            for( int i = 0; i < vmi.size(); i++ ) 
            for( int j = 0; j < vmi.size(); j++ ) 
            {
                if( i != j ) 
                    assert( vmi[i] != vmi[j] );
            }

            for( int i = 0; i < vmi.size(); i++ ) 
            {
                assert( vmi[i].absolute() == pow );
            }
            
            for( const auto& mi : vmi ) 
                LOG << mi << nl;
            LOG << nl;
            
        }
    
    }

    if(true)
    {
            
        LOG << "2. Test generator for general multiindices from 1 to ..." << nl;
        
        const std::vector<int> Ns = { -1, 0, 1, 2, 3, 4 };
        const std::vector<int> Rs = {  0, 1, 2, 3, 4, 5, 6 };
        
        
        for( const int& N : Ns ) 
        for( const int& R : Rs ) 
        {
        
            LOG << N << space << R << '$' << nl;
            
            const IndexRange bereich( 1, N );
        
            const std::vector<MultiIndex> all = generateMultiIndices( bereich, R );

            if( bereich.cardinality() == 0 and R == 0 ){
                assert( all.size() == 1 );
                assert( all[0].absolute() == 0 );
                continue;
            }

            if( bereich.cardinality() == 0 and R > 0 ){
                assert( all.size() == 0 );
                continue;
            }

            
            
            assert( all.size() == binomial_integer( bereich.cardinality() + R - 1, R ) );
            
            for( int i = 0; i < all.size(); i++ )
            for( int j = 0; j < all.size(); j++ )
                if( i != j )
                    assert( all[i] != all[j] );
            
            for( const MultiIndex& mi : all ) {
                assert( mi.absolute() == R );
            }
                
            for( const MultiIndex& mi : all )
                LOG << mi << nl;
            
            
        }
        
        LOG << "Tested" << nl;
            
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
