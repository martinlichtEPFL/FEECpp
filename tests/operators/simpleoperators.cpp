

#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../operators/simpleoperators.hpp"



int main( int argc, char *argv[] )
{
        LOG << "Unit Test for Simple operator" << nl;
        
        {
        
                LOG << "Unit Test for Scaling Component" << nl;
        
                FloatVector a(5);
                for( int i = 0; i < 5; i++ )
                        a.setentry( i, i+1 );
                
                ScalingOperator S( 5, Constants::pi );
                
                LOG << "We start with this Vector:" << nl;
                LOG << a << nl;
                
                LOG << "Scaled with PI:" << nl;
                LOG << S * a << nl;
                
                LOG << "The Scaling is " << S.getscaling() << nl;
                S.setscaling( 2.718 );
                LOG << "Now the Scaling is " << S.getscaling() << nl;
                LOG << "Accordingly:" << nl;
                LOG << S * a << nl;
                
                LOG << "Product of scaling by 11 and then 4" << nl;
                LOG << ScalingOperator(10,4) * ScalingOperator(10,11) << nl;
                LOG << ScalingOperator(10,11) * ScalingOperator(10,4) << nl;
                
                FloatVector b(10);
                for( int i = 0; i < 10; i++ )
                        b.setentry( i, i+1 );
                
                LOG << ScalingOperator(10,11) * ScalingOperator(10,4) * b << nl;
        
        }


        {
        
                LOG << "Unit Test for Diagonal Component" << nl;
        
                int dim = 10;
        
                FloatVector dia(dim);
                for( int i = 0; i < dim; i++ )
                        dia.setentry( i, i * i );
                
                DiagonalOperator D( dia );
                
                LOG << "We start with these entries:" << nl;
                LOG << dia << nl;
                
                LOG << "This is the diagonal matrix:" << nl;
                LOG << D << nl;
                
                LOG << "Pick a vector:" << nl;
                FloatVector x(dim);
                for( int i = 0; i < dim; i++ )
                        x.setentry( i, 3.01 );
                LOG << x << nl;
                LOG << "Apply the diagonal operator:" << nl;
                LOG << D * x << nl;
                
                LOG << "Now the product of diagonal operators: " << nl;
                LOG << D * D << nl;
                
                LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

        }
        
        {
        
                LOG << "Unit Test for Lambda Component" << nl;
        
                int dim = 10;
        
                FloatVector vec(dim);
                for( int i = 0; i < dim; i++ )
                        vec.setentry( i, i * i );
                
                LambdaOperator L( dim,
                        [](const FloatVector& input ) -> FloatVector
                        {
                                return 3. * input;
                        }
                );
                
                LOG << "Pick a vector:" << nl;
                FloatVector x(dim);
                for( int i = 0; i < dim; i++ )
                        x.setentry( i, i );
                LOG << x << nl;
                LOG << "Apply the Lambda operator:" << nl;
                LOG << L * x << nl;
                
                LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

        }
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

        return 0;
}
