

/**/

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.flags.hpp"
#include "../../utility/convergencetable.hpp"

#include "../../fem/global.avgsullivan.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
        
    LOG << "Unit Test: (3D) degree elevation of interpolation preserves mass" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = UnitSimplex3D();
    
    M.automatic_dirichlet_flags();
    
    M.check();
    
    
    const int r_min = 1;
    
    const int r_max = 2;
    
    const int l_min = 0;
    
    const int l_max = 3;

    const int n = M.getinnerdimension();
    
    const int number_of_samples = 3;
        
    
    Float errors[ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];

    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        for( int k = 0;     k <= n;     k++ ) 
        for( int r = r_min; r <= r_max; r++ ) 
        {
            
            LOG << "Level:       " << l_min << " <= " << l << " <= " << l_max << nl;
            
            LOG << "Polydegree:  " << r_min << " <= " << r << " <= " << r_max << nl;

            LOG << "Form degree: " << k << nl;

            LOG << "assemble mass matrices..." << nl;
            
            auto massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k, r );
            
            auto inclusion  = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), k, r );
            
            auto averaging  = FEECSullivanAveragingMatrix( M, M.getinnerdimension(), k, r, FEECAveragingMode::weighted_uniformly );
            
            auto flagmatrix = FEECSullivanFlagMatrix( M, M.getinnerdimension(), k, r );
            
            {
                auto matrix = MatrixCSR( averaging ); 

                for( int row = 0; row < matrix.getdimout(); row++ ) {
                    
                    Float sum = 0.;
                    
                    for( int i = matrix.getA()[row]; i < matrix.getA()[row+1]; i++ ) {
                        Assert( 0 <= i and i < matrix.getnumberofentries(), i );
                        sum += matrix.getV()[i];
                    }

                    if( flagmatrix.getdiagonal()[row] == 0. )
                        Assert( sum == 0., sum );
                    else 
                        Assert( is_numerically_close( sum, 1. ), sum );
                }

            }


            
            errors[k][ l-l_min ][ r-r_min ] = 0.;
            
            for( int i = 0; i < number_of_samples; i++ ){

                auto field = inclusion.createinputvector();
                field.random();

                field = flagmatrix * field;
                if ( field.norm() > 0 ) field.normalize();
                
                assert( field.isfinite() );
                
                const auto included = inclusion * field;
                
                const auto averaged = averaging * included;
                
                const auto error_eucl = ( averaged - field ).norm();
                
                const auto error_mass = ( included - inclusion * averaged ).norm(massmatrix);

                // LOG << (inclusion * averaged).norm(massmatrix) << space << included.norm(massmatrix) << space << error_eucl << space << error_mass << nl; 
                LOG << error_eucl << space << error_mass << nl; 
                
                if( error_mass >= desired_closeness )
                {
                    LOG << included - inclusion * averaged << nl;
                    exit(1);
                }

                // Last secure commit: 0ba56fe0cf92ab656d303b3fde331cb4f9b0d578
                // Still works Mar 1:  22b838bbe40d311c5675170b6f2310b296a1a8f2    
                // Still works Mar 10: 1f2d2b02a6ed96087db2d5a262b4d60a87bea166 
                // DOES NOT WORK:      a05c56458585ea03e133849b0426cf7675001923   

                Assert( error_eucl >= -desired_closeness, error_eucl ) ;
                Assert( error_mass >= -desired_closeness, error_mass ) ;
                
                Float error = error_mass;

                error = maximum( error, error_eucl );
                
                errors[k][l-l_min][r-r_min] = maximum( errors[k][l-l_min][r-r_min], error );
                
            }
            
        }
        
        if( l != l_max )
        {
            LOG << "Refinement..." << nl;
        
            M.uniformrefinement();

            LOG << "Distortion..." << nl;
        
            M.shake_interior_vertices();
            // for( auto& x : M.getcoordinates().raw() )
            // {
            //     x = sqrt(x);
            // }

        }
        
    } 

    LOG << "Convergence tables" << nl;


    std::vector<ConvergenceTable> contables( n+1 );
    
    for( int k = 0; k <= n; k++ ) 
        contables[k].table_name = "Rounding errors D3K" + std::to_string(k);
    
    for( int k = 0; k <= n; k++ ) 
    {
        for( int r = r_min; r <= r_max; r++ ) 
            contables[k] << ( "R" + std::to_string(r) );

        contables[k] << nl;     
    }

    for( int k = 0; k <= n; k++ ) 
    for( int l = l_min; l <= l_max; l++ ) 
    {
        for( int r = r_min; r <= r_max; r++ ) 
            contables[k] << errors[k][l-l_min][r-r_min];
        
        contables[k] << nl;     
    }
    
    LOG << "Check that differences are small: " << desired_closeness << nl;
    
    for( int k = 0; k <= n; k++ ) 
    {

        contables[k].lg(); 
        LOG << "-------------------" << nl;

        for( int l = l_min; l <= l_max; l++ ) 
        for( int r = r_min; r <= r_max; r++ ) 
        {
            Assert( errors[k][l-l_min][r-r_min] < desired_closeness, desired_closeness, errors[k][l-l_min][r-r_min] );
        }

    }
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
