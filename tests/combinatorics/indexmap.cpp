
#include "../../basic.hpp"
#include "../../combinatorics/indexrange.hpp"
#include "../../combinatorics/indexmap.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Index Mapping" << nl;

    const IndexRange irA(  2, 5 );
    const IndexRange irB(  3, 5 );
    const IndexRange irC( -2, 2 );
    const IndexRange irD(  0, 7 );
    const IndexRange irE(  0,-1 );

    if(true)
    {
      
        LOG << "1. Test Identity of usual interval" << nl;
        
        const IndexMap id  = identityIndexMap( irA );
        LOG << id << nl;
        
        id.check();
        
        assert( id == identityIndexMap( irA.min(), irA.max() ) );
        
        for( int a : irA ) 
        {
            assert( id.has_value_in_range( a ) );
            assert( id.preimageof( a ) == a );
            assert( id[ a ] == a );
        }
        
        assert( id.isbijective()  );
        assert( id.isinjective()  );
        assert( id.issurjective() );
        assert( id.isstrictlyascending() );
        
        assert( id.is_comparable_with( id ) );
        assert( id.is_equal_to( id ) );
        assert( not id.is_less_than( id ) );
      
    }
    
    if(true)
    {
      
        LOG << "2. Test Empty Index Map" << nl;
        
        const IndexMap leer  = identityIndexMap( irE );
        LOG << leer << nl;
        
        leer.check();
        
        assert( leer == identityIndexMap( IndexRange( 0,-1 ) ) );
        
        assert( leer.isbijective()  );
        assert( leer.isinjective()  );
        assert( leer.issurjective() );
        assert( leer.isstrictlyascending() );
        
        assert( leer.is_comparable_with( leer ) );
        assert( leer.is_equal_to( leer ) );
        assert( not leer.is_less_than( leer ) );
      
    }
    
    if(true)
    {

        LOG << "3. Test Injection and Surjection" << nl;
        
        LOG << "Injection" << nl;
        
        IndexMap inj( irB, irA, {2,3,4} );
        inj[3] = 2; inj[4] = 3; inj[5] = 4;
        assert( inj[3] == 2 && inj[4] == 3 && inj[5] == 4 );
        
        inj.check();
        
        assert( inj.isinjective() );
        assert( !inj.issurjective() );
        assert( inj.getSourceRange() == irB );
        assert( inj.getTargetRange() == irA );
        assert( inj.has_value_in_range( 4 ) );
        assert( !inj.has_value_in_range( 5 ) );
        
        LOG << "Surjection" << nl;
        
        IndexMap sur( irD, irB, {4,3,5,4,3,5,4,3} );
        sur[0] = 4; sur[1] = 3; sur[2] = 5; sur[3] = 4;
        sur[4] = 3; sur[5] = 5; sur[6] = 4; sur[7] = 3;
        assert( sur[0] == 4 && sur[1] == 3 && sur[2] == 5 && sur[3] == 4 );
        assert( sur[4] == 3 && sur[5] == 5 && sur[6] == 4 && sur[7] == 3 );

        sur.check();
        
        assert( !sur.isinjective() );
        assert( sur.issurjective() );
        assert( sur.getSourceRange() == irD );
        assert( sur.getTargetRange() == irB );
        
        LOG << "Test composition of injection with surjection" << nl;
        const IndexMap prod = inj * sur;
        LOG << inj << sur << prod << nl;
        
        prod.check();
        
        IndexMap test( irD, irA, {3,2,4,3,2,4,3,2} );
        test[0] = 3; test[1] = 2; test[2] = 4; test[3] = 3;
        test[4] = 2; test[5] = 4; test[6] = 3; test[7] = 2;
        assert( test[0] == 3 && test[1] == 2 && test[2] == 4 && test[3] == 3 );
        assert( test[4] == 2 && test[5] == 4 && test[6] == 3 && test[7] == 2 );
        
        test.check();
        
        assert( prod.is_comparable_with( test ) );
        assert( prod == test );
        
        assert( test.has_value_in_range( 3 ) );
        assert( test.preimageof( 3 ) == 0 );
        assert( test.has_value_in_range( 2 ) );
        assert( test.preimageof( 2 ) == 1 );
        assert( not test.has_value_in_range( 5 ) );
        
    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
