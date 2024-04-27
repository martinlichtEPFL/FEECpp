

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.veewedgehodge.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (2D) degree elevations commute" << nl;
        
        LOG << "Initial mesh..." << nl;
        
        auto M = UnitSquare2D();
        
        M.check();
        


        const int r_min = 0;
        
        const int r_max = 2;
        
        const int l_min = 0;
        
        const int l_max = 2;

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
        
                SparseMatrix broken_mass_matrix        = FEECBrokenMassMatrix ( M, M.getinnerdimension(),   k, r );
                SparseMatrix broken_hodge_matrix       = FEECBrokenHodgeMatrix( M, M.getinnerdimension(),   k, r );
                SparseMatrix broken_mass_hodged_matrix = FEECBrokenMassMatrix ( M, M.getinnerdimension(), 2-k, r );
                
                assert( broken_hodge_matrix.getdimin()  == broken_mass_matrix.getdimin()  );
                assert( broken_hodge_matrix.getdimout() == broken_mass_hodged_matrix.getdimout() );

                assert( broken_mass_matrix.isfinite() );
                assert( broken_hodge_matrix.isfinite() );
                
                errors[k][ l ][ r ] = 0.;
                
                for( int i = 0; i < number_of_samples; i++ ){

                    auto field = broken_hodge_matrix.createinputvector();
                    field.random();
                    field.normalize();
                    
                    assert( field.isfinite() );

                    const auto hodged_field = broken_hodge_matrix * field;
                    
                    Float mass        = field * ( broken_mass_matrix * field );
                    
                    Float hodged_mass = hodged_field * ( broken_mass_hodged_matrix * hodged_field );
                    
                    Assert( mass >= -desired_closeness, mass );
                    Assert( hodged_mass >= -desired_closeness, hodged_mass );
                    
                    const auto error_mass = absolute( mass - hodged_mass );

                    LOG << mass << space << hodged_mass << space << hodged_mass / mass << space << error_mass << nl;
                    
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
        for( int k = 0; k <= n; k++ ) {
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
