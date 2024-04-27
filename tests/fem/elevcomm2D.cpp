

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (2D) degree elevations commute" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        auto M = UnitSquare2D();
        
        M.check();
        


        const int r_min = 0;
        
        const int r_max = 5;
        
        const int l_min = 0;
        
        const int l_max = 4;

        const int n = M.getinnerdimension();
        
        const int number_of_samples = 3;
        
           
        
        Float errors[ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
            
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
            for( int k = 0; k <= n; k++ ) 
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

                LOG << "Form degree: " << space << k << nl;

                LOG << "assemble matrices..." << nl;
        
                SparseMatrix elevation_r_1 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r  , 1 );
                SparseMatrix elevation_r_2 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r+1, 1 );
                SparseMatrix elevation_r_3 = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r+2, 1 );
                SparseMatrix elevation_r_g = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k, r  , 3 );
                
                errors[k][ l ][ r ] = 0.;
                
                for( int i = 0; i < number_of_samples; i++ ){

                    auto field = elevation_r_g.createinputvector();
                    field.random();
                    field.normalize();
                    
                    assert( field.isfinite() );
                    
                    const auto path_direct   = elevation_r_g * field;
                    
                    const auto path_indirect = elevation_r_3 * elevation_r_2 * elevation_r_1 * field;
                    
                    const auto error_mass = ( path_direct - path_indirect ).norm();
                    
                    Assert( error_mass >= -desired_closeness, error_mass );
                    
                    errors[k][l-l_min][r-r_min] = maximum( errors[k][l-l_min][r-r_min], error_mass );
                    
                }
                
            }

            if( l != l_max )
            {
                LOG << "Refinement..." << nl;
            
                M.uniformrefinement();

                M.shake_interior_vertices();
            }
            
        } 
        
        
        
        LOG << "Convergence tables" << nl;
    
        std::vector<ConvergenceTable> contables( n+1 );
        
        for( int k = 0; k <= n; k++ ) 
            contables[k].table_name = "Rounding errors D2K" + std::to_string(k);
        
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
        
        
        
        for( int k = 0; k <= n; k++ ) 
        {
            contables[k].lg(); 
            LOG << "-------------------" << nl;
        }
        
        
        
        LOG << "Check that differences are small: " << desired_closeness << nl;
        
        for( int l = l_min; l <= l_max; l++ ) 
        for( int r = r_min; r <= r_max; r++ ) 
        for( int k = 0; k <= n; k++ ) 
        {
            Assert( errors[k][l-l_min][r-r_min] < desired_closeness, errors[k][l-l_min][r-r_min], desired_closeness );
        }
        
        
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
