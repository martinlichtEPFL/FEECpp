
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"
#include "../../operators/composedoperators.hpp"


int main( int argc, char *argv[] )
{
        LOG << "Unit Test for Produkt Operator Class" << nl;
        
        // if( true ) {
          
        //     LOG << "Test with two scaling matrices" << nl;
        
        //     ScalingOperator S1( 5, Constants::pi );
        //     ScalingOperator S2( 5, 4.6692 );
        //     ScalingOperator S3( 5, Constants::feigenbaum_first );
        //     ScalingOperator S4( 5, Constants::feigenbaum_second );
        //     ScalingOperator S5( 5, Constants::sirpinski );
            
        //     FloatVector x( 5 );
        //     for( int i = 0; i < 5; i++ )
        //         x[i] = i*i;
                
        //     LOG << x << nl;
            
        //     auto test1 = S5 * ( ProduktOperator( S1, S2 ) * ProduktOperator( S3, S4 ) );
        //     auto test2 = S1 * S2;
            
        //     LOG << test1 * x << nl;
        //     LOG << test2 * x << nl;
            
        // }
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

        return 0;
}
