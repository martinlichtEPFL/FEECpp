

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (3D) degree elevations commute" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        auto M = UnitCube3D();
        
        M.check();
        


        const int r_min = 0;
        
        const int r_max = 3;
        
        const int r_plus_max = 3;

        const int l_min = 0;
        
        const int l_max = 2;
        
        const int n = M.getinnerdimension();
        
        const int number_of_samples = 3;
        
           
        
        Float errors[ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        
        
            
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
            for( int k      = 0;     k <= n; k++      ) 
            for( int r      = r_min; r <= r_max;                 r++      ) 
            for( int r_plus = 0;     r_plus <= r_plus_max;       r_plus++ ) 
            {
                
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << " +" << r_plus << nl;

                LOG << "Form degree: " << space << k << nl;

                LOG << "assemble matrices..." << nl;
        
                SparseMatrix elevation = FEECBrokenElevationMatrix    ( M, M.getinnerdimension(), k, r, r_plus );
                SparseMatrix interpol  = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), k, r, r_plus );
                
                errors[k][ l ][ r ][ r_plus ] = 0.;
                
                for( int i = 0; i < number_of_samples; i++ ){

                    auto field = elevation.createinputvector();
                    field.random();
                    field.normalize();
                    
                    assert( field.isfinite() );
                    
                    const auto error_mass = ( field - interpol * elevation * field ).norm();
                    
                    Assert( error_mass >= -desired_closeness, error_mass );
                    
                    errors[k][l-l_min][r-r_min][r_plus] = maximum( errors[k][l-l_min][r-r_min][r_plus], error_mass );
                    
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
            contables[k].table_name = "Rounding errors D3K" + std::to_string(k);
        
        for( int k = 0; k <= n; k++ ) {
            for( int r = r_min; r <= r_max; r++ ) 
            for( int r_plus = 0; r_plus <= r_plus_max; r_plus++ ) 
                contables[k] << ( "R" + std::to_string(r) + "+" + std::to_string(r_plus) );

            contables[k] << nl;
        }

        for( int k = 0; k <= n; k++ ) 
        for( int l = l_min; l <= l_max; l++ ) 
        {            
            for( int r = r_min; r <= r_max; r++ ) 
            for( int r_plus = 0; r_plus <= r_plus_max; r_plus++ ) 
                contables[k] << errors[k][l-l_min][r-r_min][r_plus];
            
            contables[k] << nl;    
        }
        
        
        
        for( int k = 0; k <= n; k++ ) 
        {
            contables[k].lg(); 
            LOG << "-------------------" << nl;
        }
        
        
        
        LOG << "Check that differences are small: " << desired_closeness << nl;
        
        for( int k      =     0; k      <= n; k++      ) 
        for( int l      = l_min; l      <=                 l_max; l++      ) 
        for( int r      = r_min; r      <=                 r_max; r++      ) 
        for( int r_plus =     0; r_plus <=            r_plus_max; r_plus++ ) 
        {
            Assert( errors[k][l-l_min][r-r_min][r_plus] < desired_closeness, errors[k][l-l_min][r-r_min][r_plus], desired_closeness );
        }
        
        
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
