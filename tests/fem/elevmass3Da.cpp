

/**/

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (3D) degree elevation of interpolation preserves mass" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.check();
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
        std::vector<Float>                                          experiments_scalar_value;
        
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                     assert( vec.getdimension() == 3 );
//                     return FloatVector({ 1. });
//                 }
//         );
// 
//         experiments_scalar_value.push_back( 1. );
//         
//         
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                     assert( vec.getdimension() == 3 );
//                     return FloatVector({ vec.sum() < 1. ? 1. : 0. });
//                 }
//         );
// 
//         experiments_scalar_value.push_back( 1./6. );
        
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ std::exp( vec[0] + vec[1] + vec[2] ) });
            }
        );

        experiments_scalar_value.push_back( 3.19452804946532511361521373028750 * 3.19452804946532511361521373028750 * 3.19452804946532511361521373028750 );


        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_field;
        std::vector<Float>                                          experiments_vector_value;
        
//         experiments_vector_field.push_back( 
//         [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 3 );
//                 return FloatVector({ 1.,0.,0. });
//             }
//         );
// 
//         experiments_vector_value.push_back( 1. );
//         
//         
//         experiments_vector_field.push_back( 
//         [](const FloatVector& vec) -> FloatVector{
//                 assert( vec.getdimension() == 3 );
//                 return FloatVector({ 1.,-3.,-2. });
//             }
//         );
//         
//         experiments_vector_value.push_back( 14. );
        
        
        experiments_vector_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    Float x = vec[0]; Float y = vec[1]; Float z = vec[2]; 
                    return FloatVector({ std::sin( x*y ), std::cos( z ), std::exp(x+y+z) });
                }
            );

        experiments_vector_value.push_back( 33.42616007376754121867328484399462245509699700639766495311 );
        
        
        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_pseudo_field;
        std::vector<Float>                                          experiments_pseudo_value;

        experiments_pseudo_field = experiments_vector_field;
        experiments_pseudo_value = experiments_vector_value;
        
        
        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
        std::vector<Float>                                          experiments_volume_value;

        experiments_volume_field = experiments_scalar_field;
        experiments_volume_value = experiments_scalar_value;

        
        
        const int r_min = 0;
        
        const int r_max = 2;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        const int r_plus_max = 3;
        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        Float errors_vector[ experiments_vector_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        Float errors_pseudo[ experiments_pseudo_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
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
                
                SparseMatrix massmatrix_pseudo      = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
                
                SparseMatrix massmatrix_volume      = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r );
                
                SparseMatrix massmatrix_scalar_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r + r_plus );
                
                SparseMatrix massmatrix_vector_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + r_plus );
                
                SparseMatrix massmatrix_pseudo_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r + r_plus );
                
                SparseMatrix massmatrix_volume_plus = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r + r_plus );
                
                LOG << "assemble degree elevation matrices..." << nl;

                SparseMatrix elevation_scalar       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus );

                SparseMatrix elevation_vector       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r, r_plus );

                SparseMatrix elevation_pseudo       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 2, r, r_plus );

                SparseMatrix elevation_volume       = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 3, r, r_plus );

                
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
                
                for( int i = 0; i < experiments_pseudo_field.size(); i++ ){

                    const auto& pseudofield = experiments_pseudo_field[i];
        
                    FloatVector interpol      = Interpolation( M, M.getinnerdimension(), 2, r, pseudofield );

                    FloatVector interpol_elev = elevation_pseudo * interpol;
                    
                    Float mass      = interpol * ( massmatrix_pseudo * interpol );

                    Float mass_elev = interpol_elev * ( massmatrix_pseudo_plus * interpol_elev );

                    Assert( mass >= -desired_closeness, mass );
                    Assert( mass_elev >= -desired_closeness, mass_elev);
                    
                    Float error_mass = mass - mass_elev;
                    
                    errors_pseudo[i][l-l_min][r-r_min][r_plus] = error_mass;
                    
                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++ ){

                    const auto& volumefield = experiments_volume_field[i];
        
                    FloatVector interpol      = Interpolation( M, M.getinnerdimension(), 3, r, volumefield );

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
        ConvergenceTable contable_pseudo[ experiments_vector_field.size() ];
        ConvergenceTable contable_volume[ experiments_volume_field.size() ];
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i].table_name = "Rounding errors scalar E" + std::to_string(i);
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i].table_name = "Rounding errors vector E" + std::to_string(i);
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                contable_pseudo[i].table_name = "Rounding errors pseudo E" + std::to_string(i);
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i].table_name = "Rounding errors volume E" + std::to_string(i);

            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i] << printf_into_string("R%d+%d", r-r_min, r_plus_max );
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i] << printf_into_string("R%d+%d", r-r_min, r_plus_max );
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                contable_pseudo[i] << printf_into_string("R%d+%d", r-r_min, r_plus_max );
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i] << printf_into_string("R%d+%d", r-r_min, r_plus_max );

        }
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
        for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i] << nl; 
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
        
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                    contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min][r_plus_max];
            
                for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                    contable_vector[i] << errors_vector[i][l-l_min][r-r_min][r_plus_max];
            
                for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                    contable_pseudo[i] << errors_pseudo[i][l-l_min][r-r_min][r_plus_max];
            
                for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                    contable_volume[i] << errors_volume[i][l-l_min][r-r_min][r_plus_max];
            
            }
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
            for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i] << nl; 
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
        
        LOG << "Convergence tables: pseudos" << nl;
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
        {
            contable_pseudo[i].lg(); 
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
            
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                Assert( errors_pseudo[i][l-l_min][r-r_min][r_plus] < desired_closeness, errors_pseudo[i][l-l_min][r-r_min][r_plus], desired_closeness );
            
            for( int i = 0; i < experiments_volume_field.size(); i++ )
                Assert( errors_volume[i][l-l_min][r-r_min][r_plus] < desired_closeness, errors_volume[i][l-l_min][r-r_min][r_plus], desired_closeness );
        }
            
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
