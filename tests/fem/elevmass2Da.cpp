

/**/

#include "../../basic.hpp"
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
        
        LOG << "Unit Test: (2D) degree elevation of interpolation preserves mass" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ std::exp( vec[0] ) });
            }
        );

        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                auto x = vec[0];
                return FloatVector({ 
                    ( x*x - 1. ) / ( x*x + 1. ) 
                });
            }
        );

        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_field;
        
        experiments_vector_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ 
                    std::exp( vec[0] ), 
                    std::sin( vec[1] )
                });
            }
        );

        experiments_vector_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                auto x = vec[0];
                auto y = vec[1];
                return FloatVector({ 
                    ( x*x - 1. ) / ( y*y       + 1. ),
                    ( x*y + 1. ) / ( y*y + x*x + 1. )                
                });
            }
        );

        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
        
        experiments_volume_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ 
                    std::exp( vec[0] )
                });
            }
        );
        
        experiments_volume_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                auto x = vec[0];
                return FloatVector({ 
                    ( x*x - 1. ) / ( x*x + 1. )
                });
            }
        );

        
        
        const int r_min = 0;
        
        const int r_max = 3;
        
        const int l_min = 0;
        
        const int l_max = 3;
        
        const int r_plus_max = 3;
        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        Float errors_vector[ experiments_vector_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
            for( int r      = r_min; r      <=      r_max; r++      ) 
            for( int r_plus =     0; r_plus <= r_plus_max; r_plus++ ) 
            {
                
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;
                LOG << "Adding: 0 <= " << r_plus << " <= " << r_plus_max << nl;

                LOG << "assemble mass matrices..." << nl;

                SparseMatrix massmatrix_scalar      = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                
                SparseMatrix massmatrix_vector      = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                
                SparseMatrix massmatrix_volume      = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
                
                SparseMatrix massmatrix_scalar_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r + r_plus );
                
                SparseMatrix massmatrix_vector_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + r_plus );
                
                SparseMatrix massmatrix_volume_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r + r_plus );
                
                LOG << "assemble degree elevation matrices..." << nl;

                SparseMatrix elevation_scalar       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus );

                SparseMatrix elevation_vector       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_plus );

                SparseMatrix elevation_volume       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r, r_plus );

                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ){

                    const auto& scalarfield = experiments_scalar_field[i];
        
                    FloatVector interpol      = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                    FloatVector interpol_elev = elevation_scalar * interpol;
                    
                    Float mass      = interpol * ( massmatrix_scalar * interpol );

                    Float mass_elev = interpol_elev * ( massmatrix_scalar_plus * interpol_elev );

                    Assert( mass >= -desired_closeness, mass );
                    Assert( mass_elev >= -desired_closeness, mass_elev);
                    
                    Float error_mass = mass - mass_elev;
                    
                    errors_scalar[i][l-l_min][r-r_min][r_plus] = error_mass;
                    
                }
                
                for( int i = 0; i < experiments_vector_field.size(); i++ ){

                    const auto& vectorfield = experiments_vector_field[i];
        
                    FloatVector interpol      = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );

                    FloatVector interpol_elev = elevation_vector * interpol;
                    
                    Float mass      = interpol * ( massmatrix_vector * interpol );

                    Float mass_elev = interpol_elev * ( massmatrix_vector_plus * interpol_elev );

                    Assert( mass >= -desired_closeness, mass );
                    Assert( mass_elev >= -desired_closeness, mass_elev);
                    
                    Float error_mass = mass - mass_elev;
                    
                    errors_vector[i][l-l_min][r-r_min][r_plus] = error_mass;
                    
                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++ ){

                    const auto& volumefield = experiments_volume_field[i];
        
                    FloatVector interpol      = Interpolation( M, M.getinnerdimension(), 2, r, volumefield );

                    FloatVector interpol_elev = elevation_volume * interpol;
                    
                    Float mass      = interpol * ( massmatrix_volume * interpol );

                    Float mass_elev = interpol_elev * ( massmatrix_volume_plus * interpol_elev );

                    Assert( mass >= -desired_closeness, mass );
                    Assert( mass_elev >= -desired_closeness, mass_elev);
                    
                    Float error_mass = mass - mass_elev;
                    
                    errors_volume[i][l-l_min][r-r_min][r_plus] = error_mass;
                    
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
    
        ConvergenceTable contable_scalar[ experiments_scalar_field.size() ];
        ConvergenceTable contable_vector[ experiments_vector_field.size() ];
        ConvergenceTable contable_volume[ experiments_volume_field.size() ];
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i].table_name = "Rounding errors scalar E" + std::to_string(i);
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i].table_name = "Rounding errors vector E" + std::to_string(i);
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i].table_name = "Rounding errors volume E" + std::to_string(i);

            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i] << printf_into_string("R%d+%d", r-r_min, r_plus_max );
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i] << printf_into_string("R%d+%d", r-r_min, r_plus_max );
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i] << printf_into_string("R%d+%d", r-r_min, r_plus_max );

        }
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
        for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                    contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min][r_plus_max];
            
                for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                    contable_vector[i] << errors_vector[i][l-l_min][r-r_min][r_plus_max];
            
                for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                    contable_volume[i] << errors_volume[i][l-l_min][r-r_min][r_plus_max];
            
            }
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
            for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
            for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
            
        }
        
        
        
        LOG << "Convergence tables: scalars" << nl;
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
        {
            contable_scalar[i].lg(); 
            LOG << "-------------------" << nl;
        }
        
        LOG << "Convergence tables: vectors" << nl;
        for( int i = 0; i < experiments_vector_field.size(); i++ ) 
        {
            contable_vector[i].lg(); 
            LOG << "-------------------" << nl;
        }
        
        LOG << "Convergence tables: volumes" << nl;
        for( int i = 0; i < experiments_volume_field.size(); i++ )
        {
            contable_volume[i].lg(); 
            LOG << "-------------------" << nl;
        }
        
        
        
        
        
        LOG << "Check that differences are small: " << desired_closeness << nl;
        
        for( int l      = l_min; l      <=      l_max; l++      ) 
        for( int r      = r_min; r      <=      r_max; r++      ) 
        for( int r_plus =     0; r_plus <= r_plus_max; r_plus++ ) 
        {
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                Assert( errors_scalar[i][l-l_min][r-r_min][r_plus] < desired_closeness, errors_scalar[i][l-l_min][r-r_min][r_plus], desired_closeness );
            
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                Assert( errors_vector[i][l-l_min][r-r_min][r_plus] < desired_closeness, errors_vector[i][l-l_min][r-r_min][r_plus], desired_closeness );
            
            for( int i = 0; i < experiments_volume_field.size(); i++ )
                Assert( errors_volume[i][l-l_min][r-r_min][r_plus] < desired_closeness, errors_volume[i][l-l_min][r-r_min][r_plus], desired_closeness );
        }
            
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
