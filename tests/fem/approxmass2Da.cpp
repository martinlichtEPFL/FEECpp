

/**/

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (2D) masses are correctly approximated: precomputed mass" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial2D M = StandardSquare2D();
        
        M.check();
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
        std::vector<Float>                                          experiments_scalar_value;
        
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 Assert( vec.getdimension() == 2 );
//                 return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? 1. : Float(0.) });
//             }
//         );
// 
//         experiments_scalar_value.push_back( 1. );
//         
//         
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 Assert( vec.getdimension() == 2 );
//                 return FloatVector({ 1. });
//             }
//         );
// 
//         experiments_scalar_value.push_back( 4. );
//         
//         
//         experiments_scalar_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 Assert( vec.getdimension() == 2 );
//                 return FloatVector({ ( vec[0] * vec[1] > 0 ) ? 1. : Float(0.) });
//             }
//         );
// 
//         experiments_scalar_value.push_back( 2. );
        
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                Assert( vec.getdimension() == 2 );
                return FloatVector({ std::exp( vec[0] + vec[1] ) });
            }
        );

        experiments_scalar_value.push_back( 3.62686040784701876 * 3.62686040784701876 );

        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                Assert( vec.getdimension() == 2 );
                return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? std::exp( vec[0] ) : Float(0.) });
            }
        );

        experiments_scalar_value.push_back(3.194528049465325113615213730287503906590157785275923662043);
        
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                Assert( vec.getdimension() == 2 );
                return FloatVector({ ( vec[0] * vec[1] > 0 ) ? std::exp( vec[0] ) : Float(0.) });
            }
        );

        experiments_scalar_value.push_back(3.626860407847018767668213982801261704886342012321135721309);
        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_field;
        std::vector<Float>                                          experiments_vector_value;
        
//         experiments_vector_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 Assert( vec.getdimension() == 2 );
//                 // return FloatVector({ ( vec[0]*vec[1] > 0 ) ? 1. : 0., ( vec[0]*vec[1] > 0 ) ? 10. : Float(0.) });
//                 return FloatVector({ 1., Float(0.) });
//             }
//         );
// 
//         experiments_vector_value.push_back( 4. );
        
        
        experiments_vector_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                Assert( vec.getdimension() == 2 );
                Float x = vec[0]; Float y = vec[1];
                return FloatVector({ std::exp(x-y), std::exp(x-y) });
            }
        );

        // int_(-1)^(1) int_(-1)^(1) exp(x-y)^2 + exp(x-y)^2 dx dy
        experiments_vector_value.push_back( 26.308232836016486629201989612067059822501324553083772160298096942 );
        

        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
        std::vector<Float>                                          experiments_volume_value;

//         experiments_volume_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 Assert( vec.getdimension() == 2 );
//                 return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? 1. : Float(0.) });
//             }
//         );
// 
//         experiments_volume_value.push_back( 1. );
//         
//         
//         experiments_volume_field.push_back( 
//             [](const FloatVector& vec) -> FloatVector{
//                 Assert( vec.getdimension() == 2 );
//                 return FloatVector({ ( vec[0] * vec[1] > 0 ) ? 1. : Float(0.) });
//             }
//         );
// 
//         experiments_volume_value.push_back( 2. );
        
        
        experiments_volume_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                Assert( vec.getdimension() == 2 );
                return FloatVector({ ( vec[0] > 0 and vec[1] > 0 ) ? std::exp( vec[0] ) : Float(0.) });
            }
        );

        experiments_volume_value.push_back(3.194528049465325113615213730287503906590157785275923662043);
        
        
        experiments_volume_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                Assert( vec.getdimension() == 2 );
                return FloatVector({ ( vec[0] * vec[1] > 0 ) ? std::exp( vec[0] ) : Float(0.) });
            }
        );

        experiments_volume_value.push_back(3.626860407847018767668213982801261704886342012321135721309);




        
        
        const int r_min = 0;
        
        const int r_max = 5;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_vector[ experiments_vector_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
                    
            for( int r = r_min; r <= r_max; r++ ) 
            {
        
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

                LOG << "assemble mass matrices..." << nl;
                
                SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                
                SparseMatrix massmatrix_vector = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                
                SparseMatrix massmatrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
                
                Assert( massmatrix_scalar.isfinite() );
                Assert( massmatrix_vector.isfinite() );
                Assert( massmatrix_volume.isfinite() );
                
                LOG << "experiments..." << nl;
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ){

                    const auto& scalarfield = experiments_scalar_field[i];
                    const auto& should_be   = experiments_scalar_value[i];

                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                    Float mass = interpol * ( massmatrix_scalar * interpol );
                    
                    Assert( mass >= -desired_closeness, mass );
                    
                    errors_scalar[i][l-l_min][r-r_min] = std::sqrt( std::abs( mass - should_be ) );
                    
                }
                
                for( int i = 0; i < experiments_vector_field.size(); i++ ){

                    const auto& vectorfield = experiments_vector_field[i];
                    const auto& should_be   = experiments_vector_value[i];

                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );

                    Float mass = interpol * ( massmatrix_vector * interpol );
                    
                    Assert( mass >= -desired_closeness, mass );
                    
                    errors_vector[i][l][r] = std::sqrt( std::abs( mass - should_be ) );
                    
                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++ ){

                    const auto& volumefield = experiments_volume_field[i];
                    const auto& should_be   = experiments_volume_value[i];

                    FloatVector interpol = Interpolation( M, M.getinnerdimension(), 2, r, volumefield );

                    Float mass = interpol * ( massmatrix_volume * interpol );
                    
                    Assert( mass >= -desired_closeness, mass );
                    
                    errors_volume[i][l][r] = std::sqrt( std::abs( mass - should_be ) );
                    
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
                contable_scalar[i].table_name = "Numerical errors scalar E" + std::to_string(i);
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i].table_name = "Numerical errors vector E" + std::to_string(i);
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i].table_name = "Numerical errors volume E" + std::to_string(i);

            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i] << printf_into_string("R%d", r-r_min );
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i] << printf_into_string("R%d", r-r_min );
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i] << printf_into_string("R%d", r-r_min );

        }
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
        for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
        
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                    contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min];
            
                for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                    contable_vector[i] << errors_vector[i][l-l_min][r-r_min];
            
                for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                    contable_volume[i] << errors_volume[i][l-l_min][r-r_min];
            
            }
            
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
            for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
            for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
            
        }
            
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i].lg(); 
        LOG << "-------------------" << nl;
        for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i].lg(); 
        LOG << "-------------------" << nl;
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i].lg(); 
        
        
        
        
        
        LOG << "Check that differences are small: " << desired_closeness << nl;
        
        for( int l      = l_min; l      <=      l_max; l++      ) 
        for( int r      = r_min; r      <=      r_max; r++      ) 
        {
            if( r < r_max or l < 3 ) 
                continue;
            
            continue; // TODO: find a meaningful test here 

            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                Assert( errors_scalar[i][l-l_min][r-r_min] < desired_closeness, errors_scalar[i][l-l_min][r-r_min], desired_closeness );
            
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                Assert( errors_vector[i][l-l_min][r-r_min] < desired_closeness, errors_vector[i][l-l_min][r-r_min], desired_closeness );
            
            for( int i = 0; i < experiments_volume_field.size(); i++ )
                Assert( errors_volume[i][l-l_min][r-r_min] < desired_closeness, errors_volume[i][l-l_min][r-r_min], desired_closeness );
        }
        
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
