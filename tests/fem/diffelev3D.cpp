

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../utility/convergencetable.hpp"

#include "../../fem/global.unphysical.hpp"



using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (3D) degree elevation commutes with exterior derivative" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = StandardCube3D();
        
        M.check();
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> fields;
        
        fields.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ std::exp( vec[0] + vec[1] - vec[2] ) });
            }
        );

        fields.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                auto ret = FloatVector({ 
                                std::exp( vec[0] + vec[2] ), 
                                std::sin( -5.*vec[0] -vec[1] ),
                                std::exp( vec[1] * vec[2] )
                });
                assert( ret.getdimension() == 3 );
                return ret;
            }
        );

        fields.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                auto ret = FloatVector({ 
                                std::cos( -5.*vec[0]*vec[0] ), 
                               -std::exp( 5.*vec[2] / ( 1. + vec[0] * vec[0] ) ), 
                                std::cos( 2.*vec[1]-vec[0] ) 
                                
                });
                assert( ret.getdimension() == 3 );
                return ret;
            }
        );
        
        

        const int r_min = 1;
        
        const int r_max = 1;
        
        const int l_min = 0;
        
        const int l_max = 1;
        
        const int r_plus_max = 1;
         
        Float errors[ M.getinnerdimension() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ][ r_plus_max + 1 ];
        
        
            
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
            for( int k      =     0; k      <  M.getinnerdimension(); k++      ) 
            for( int r      = r_min; r      <=                 r_max; r++      ) 
            for( int r_plus =     0; r_plus <=            r_plus_max; r_plus++ ) 
            {
                
                LOG << "...assemble matrices: l=" << l << " k=" << k << " r=" << r << " rplus=" << r_plus << nl;
        
                SparseMatrix lower_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r          );

                SparseMatrix upper_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), k, r + r_plus );

                SparseMatrix diyi_elevation = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k  , r  , r_plus );
                
                SparseMatrix dier_elevation = FEECBrokenElevationMatrix( M, M.getinnerdimension(), k+1, r-1, r_plus );
                
                SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k+1, r + r_plus - 1 );
                
                SparseMatrix origin_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), k, r );
                
                
                
                FloatVector interpol_function = lower_diffmatrix.createinputvector();
                interpol_function.random();
                interpol_function.normalize(origin_massmatrix);

                auto path1 = dier_elevation * lower_diffmatrix * interpol_function;

                auto path2 = upper_diffmatrix * diyi_elevation * interpol_function;


                // SparseMatrix canon = FEECRandomizeBroken( M, M.getinnerdimension(), k+1, r + r_plus - 1, notanumber );
                SparseMatrix canon = FEECCanonicalizeBroken( M, M.getinnerdimension(), k+1, r + r_plus - 1 );
                
                auto commutator_error = path2 - path1;

                // LOG << massmatrix * ( canon * ( path1 - path2 ) - ( path1 - path2 ) ) << nl;

                // const auto csr_product = MatrixCSR(canon.getTranspose()) & MatrixCSR(massmatrix) & MatrixCSR(canon);
                
                // Float commutator_error_mass = commutator_error * ( massmatrix * commutator_error );
                // Float commutator_error_mass = ( canon * commutator_error ) * ( massmatrix * commutator_error );
                // Float commutator_error_mass = norm_sq_of_vector( massmatrix, canon * commutator_error );
                Float commutator_error_mass = norm_sq_of_vector( massmatrix, commutator_error );
                // auto foo = canon * commutator_error; Float commutator_error_mass = foo * ( massmatrix * foo );
                // Float commutator_error_mass = commutator_error * ( csr_product * commutator_error );

                assert( std::isfinite( commutator_error_mass ) );
                Assert( commutator_error_mass >= -desired_closeness, commutator_error_mass );
                                
                errors[k][l-l_min][r-r_min][r_plus] = std::sqrt( std::fabs( commutator_error_mass ) );
            
                
                
            }

            if( l != l_max )
            {
                LOG << "Refinement..." << nl;
            
                M.uniformrefinement();

                LOG << M.getcoordinates().text() << nl;

                M.shake_interior_vertices();

                LOG << M.getcoordinates().text() << nl;
            }
            
            

        } 
        
        
        LOG << "Convergence tables" << nl;
    
        ConvergenceTable contables[ M.getinnerdimension() ];
        
        for( int k = 0; k < M.getinnerdimension(); k++ ) 
            contables[k].table_name = "Rounding errors D3K" + std::to_string(k);
        for( int k = 0; k < M.getinnerdimension(); k++ ) {
            for( int r = r_min; r <= r_max; r++ ) 
                contables[k] << printf_into_string("R%d+%d", r-r_min, r_plus_max );;
            contables[k] << nl;
        }

        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < M.getinnerdimension(); i++ ) 
                    contables[i] << errors[i][l-l_min][r-r_min][r_plus_max];
                        
            }
            
            for( int i = 0; i < M.getinnerdimension(); i++ ) 
                contables[i] << nl; 
            
        }
        
        
        
        for( int i = 0; i < M.getinnerdimension(); i++ ) 
        {
            contables[i].lg(); 
            LOG << "-------------------" << nl;
        }
                
        
        
        
        
        LOG << "Check that differences are small: " << desired_closeness << nl;
        
        for( int l      = l_min; l      <=           l_max; l++      ) 
        for( int r      = r_min; r      <=           r_max; r++      ) 
        for( int r_plus =     0; r_plus <=      r_plus_max; r_plus++ ) 
        for( int i      =     0; i < M.getinnerdimension(); i++      ) 
        {
            Assert( errors[i][l-l_min][r-r_min][r_plus] < desired_closeness, errors[i][l-l_min][r-r_min][r_plus], desired_closeness );
        }
            
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
