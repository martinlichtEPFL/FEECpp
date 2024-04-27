

#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"



int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Index Ranges" << nl;
    
    if( true ) 
    {
        
        LOG << "1. Test empty index ranges" << nl;
        
        IndexRange irE1(  3,  2 );
        IndexRange irE2(  5,  1 );
        IndexRange irE3( -2, -3 );
        
        assert( irE1.isempty() );
        assert( irE2.isempty() );
        assert( irE3.isempty() );
        
        assert( irE1.cardinality() == 0 );
        assert( irE2.cardinality() == 0 );
        assert( irE3.cardinality() == 0 );
        
        int counter = 0;
        for( int i : irE1 ) { counter++; assert( irE1.min() <= i && i <= irE1.max() ); unreachable(); }
        for( int i : irE2 ) { counter++; assert( irE2.min() <= i && i <= irE2.max() ); unreachable(); }
        for( int i : irE3 ) { counter++; assert( irE3.min() <= i && i <= irE3.max() ); unreachable(); }
        assert( counter == 0 );
        
        assert( irE1 == irE2 );
        assert( irE1 == irE3 );
        assert( irE2 == irE3 );
        
    } 
        
    if( true ) 
    {
        
        LOG << "2. Test non-empty index ranges" << nl;
        
        IndexRange irA( 3, 7 );
        IndexRange irB( 5, 5 );
        IndexRange irC(-2,-2 );
        IndexRange irD( 4, 9 );
        
        assert( not irA.isempty() );
        assert( not irB.isempty() );
        assert( not irC.isempty() );
        assert( not irD.isempty() );
        
        assert( irA.cardinality() == 5 );
        assert( irB.cardinality() == 1 );
        assert( irC.cardinality() == 1 );
        assert( irD.cardinality() == 6 );
        
        assert( irA != irB );
        assert( irA != irC );
        assert( irA != irD );
        assert( irB != irC );
        assert( irB != irD );
        assert( irC != irD );
        
        assert( !irA.contains(2) );
        assert( irA.contains(3) );
        assert( irA.contains(4) );
        assert( irA.contains(5) );
        assert( irA.contains(6) );
        assert( irA.contains(7) );
        assert( !irA.contains(8) );
        
        assert( irB.contains(5) );
        assert( !irB.contains(4) );
        assert( !irB.contains(6) );
        
        LOG << "Test indexing in non-empty index ranges" << nl;
        
        assert( irB.element2position(5) == 0 );
        assert( irB.position2element(0) == 5 );
        assert( irA.element2position(3) == 0 );
        assert( irA.element2position(5) == 2 );
        assert( irA.element2position(7) == 4 );
        assert( irA.position2element(0) == 3 );
        assert( irA.position2element(2) == 5 );
        assert( irA.position2element(4) == 7 );
        
    }
    
    if( true )
    {
        
        LOG << "3. Test Index Range iterators over some normal intervals" << nl;
        
        std::vector<IndexRange> irs = {
            IndexRange(  2, 5 ),
            IndexRange(  7, 7 ),
            IndexRange( -3, 3 ),
            IndexRange( -6, 0 )
        };

        for( auto ir : irs ){

            assert( ! ir.isempty() );
        
            LOG << "For each loop " << ir.min() << space << ir.max() << nl << tab;
            for( int i : ir ) {
                LOG << i << space;
                assert( ir.min() <= i && i <= ir.max() );
            }
                
            LOG << nl;

            LOG << "Classical For loop " << ir.min() << space << ir.max() << nl << tab;
            for( IndexRange::ConstIterator iri = ir.begin(); iri != ir.end(); ++iri ) {
                LOG << *iri << space;
                assert( ir.min() <= *iri && *iri <= ir.max() );
            }
                
            LOG << nl;
            
            LOG << "While Loop " << ir.min() << space << ir.max() << nl << tab;
            IndexRange::ConstIterator iri = ir.begin();
            while( iri != ir.end() ) {
                int i = *(iri++);
                assert( ir.min() <= i && i <= ir.max() );
                LOG << i << space;
            }

            LOG << nl;
            
        }
        
    }
        
    if( true )
    {
    
        LOG << "4. Test Index Range iterators over some extreme cases" << nl;
        
        std::vector<IndexRange> irs = {
            IndexRange( std::numeric_limits<int>::max()   , std::numeric_limits<int>::max()    ),
            IndexRange( std::numeric_limits<int>::max()- 1, std::numeric_limits<int>::max()    ),    
            IndexRange( std::numeric_limits<int>::max()-10, std::numeric_limits<int>::max()    ),
            IndexRange( std::numeric_limits<int>::min(),    std::numeric_limits<int>::min()    ),
            IndexRange( std::numeric_limits<int>::min(),    std::numeric_limits<int>::min()+ 1 ),
            IndexRange( std::numeric_limits<int>::min(),    std::numeric_limits<int>::min()+10 ) 
        };

        for( auto ir : irs ){

            assert( ! ir.isempty() );
        
            LOG << "For each loop " << ir.min() << space << ir.max() << nl << tab;
            for( int i : ir ) {
                LOG << i << space;
                assert( ir.min() <= i && i <= ir.max() );
            }
                
            LOG << nl;

            LOG << "Classical For loop " << ir.min() << space << ir.max() << nl << tab;
            for( IndexRange::ConstIterator iri = ir.begin(); iri != ir.end(); ++iri ) {
                LOG << *iri << space;
                assert( ir.min() <= *iri && *iri <= ir.max() );
            }
                
            LOG << nl;
            
            LOG << "While Loop " << ir.min() << space << ir.max() << nl << tab;
            IndexRange::ConstIterator iri = ir.begin();
            while( iri != ir.end() ) {
                int i = *(iri++);
                assert( ir.min() <= i && i <= ir.max() );
                LOG << i << space;
            }

            LOG << nl;
            
        }
        
    }
        
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
