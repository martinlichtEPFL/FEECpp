

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../utility/math.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
    LOG << "Unit Test: 3D Poisson Problem" << nl;
        
        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial3D M = FicheraCorner3D();

            M.automatic_dirichlet_flags();
            
            M.check();
            
            M.check_dirichlet_flags();
            
            LOG << "Prepare scalar fields for testing..." << nl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 3 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    Float r; 
                    Float theta;
                    cartesian_to_polar_coordinates2D( vec[0], vec[1], r, theta );
                    const Float pi = Constants::pi; 

                    Float u_p = (-1)*r*r/(6*pi) 
                                * 
                                ( 3 * pi / 2 + 2 * log( r ) * sin( 2 * theta ) + ( 2 * theta - 3 / 2 * pi ) * cos(2*theta) );

                    const int j = 1; // only odd indices count ... 
                    const Float alpha_j = 0.40192487;
                    Float u = alpha_j * power_numerical( r, 2 * j / 3. ) * sin( 2./3. * j * theta );

                    return FloatVector({ u + u_p });
                };
        
        
            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    return FloatVector({
                        1.0
                     });
                };
            

            

            

            const int min_l = 1; 
            const int max_l = 3;
            

            ConvergenceTable contable("Mass error");
            
            contable << "u_error" << "du_error" << "residual" << "time" << nl;

            
            const int r = 1;

            const int r_plus = 3;
            
            assert( 0 <= min_l and min_l <= max_l );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << nl;
                LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                LOG << "...assemble matrices" << nl;
        
                SparseMatrix aug_scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r + r_plus     );
                
                SparseMatrix aug_vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r + r_plus - 1 );
                
                SparseMatrix aug_diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r + r_plus );

                SparseMatrix aug_incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r + r_plus );
                SparseMatrix     incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r          );

                SparseMatrix elevmatrix = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, r_plus ); 

                SparseMatrix aug_diffmatrix_t = aug_diffmatrix.getTranspose();
                SparseMatrix aug_incmatrix_t  = aug_incmatrix.getTranspose();
                
                SparseMatrix incmatrix_t  = incmatrix.getTranspose();
                SparseMatrix elevmatrix_t = elevmatrix.getTranspose(); 
                
                
                LOG << "...compose matrices" << nl;
        
                auto opr  = aug_diffmatrix;
                auto opl  = opr.getTranspose(); 
                auto opm = opl & ( aug_vector_massmatrix & opr );

                auto aug_stiffness = aug_incmatrix_t & opm & aug_incmatrix;
                
                auto stiffness = incmatrix_t & elevmatrix_t & opm & elevmatrix & incmatrix;
                
                aug_stiffness.sortentries();
                stiffness.sortentries();

                auto     stiffness_csr = MatrixCSR(     stiffness );
                auto aug_stiffness_csr = MatrixCSR( aug_stiffness );
                
                
                LOG << "...interpolate" << nl;

                const auto& function_sol  = experiment_sol;
                const auto& function_rhs  = experiment_rhs;
                
                FloatVector     interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,        function_sol  );
                FloatVector     interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,        function_rhs  );
                FloatVector aug_interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r+r_plus, function_rhs  );

                FloatVector     rhs =     incmatrix_t * elevmatrix_t * aug_scalar_massmatrix * elevmatrix * interpol_rhs;
                FloatVector aug_rhs = aug_incmatrix_t * aug_scalar_massmatrix * aug_interpol_rhs;

                FloatVector     sol(     incmatrix.getdimin(), 0. );
                FloatVector aug_sol( aug_incmatrix.getdimin(), 0. );
                
                timestamp start = timestampnow();

                {
                    LOG << "...iterative solver 1" << nl;
                
                    sol.zero();
                    
                    auto diagonal = stiffness_csr.getDiagonal();
                    FloatVector residual( rhs );

                    ConjugateGradientSolverCSR_SSOR( 
                        sol.getdimension(), 
                        sol.raw(), 
                        rhs.raw(), 
                        stiffness_csr.getA(), stiffness_csr.getC(), stiffness_csr.getV(),
                        residual.raw(),
                        desired_precision,
                        0,
                        diagonal.raw(),
                        1.0
                    );

                }

                {
                    LOG << "...iterative solver augmented" << nl;
                
                    aug_sol.zero();
                    
                    auto aug_diagonal = aug_stiffness_csr.getDiagonal();
                    FloatVector aug_residual( aug_rhs );

                    ConjugateGradientSolverCSR_SSOR( 
                        aug_sol.getdimension(), 
                        aug_sol.raw(), 
                        aug_rhs.raw(), 
                        aug_stiffness_csr.getA(), aug_stiffness_csr.getC(), aug_stiffness_csr.getV(),
                        aug_residual.raw(),
                        desired_precision,
                        0,
                        aug_diagonal.raw(),
                        1.0
                    );

                }

                timestamp end = timestampnow();
                LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                LOG << "...compute error and residual" << nl;

                FloatVector error     = aug_incmatrix * aug_sol - elevmatrix * incmatrix * sol;
                FloatVector graderror = aug_diffmatrix * ( aug_incmatrix * aug_sol - elevmatrix * incmatrix * sol );
                Float errornorm       = std::sqrt( error * ( aug_scalar_massmatrix * error ) );
                Float graderrornorm   = std::sqrt( graderror * ( aug_vector_massmatrix * graderror ) );
                Float residualnorm  = ( rhs - stiffness * sol ).norm();

                LOG << "error:     " << errornorm    << nl;
                LOG << "graderror: " << graderrornorm << nl;
                LOG << "residual:  " << residualnorm << nl;
                LOG << "time:      " << Float( end - start ) << nl;
                        
                contable << errornorm;
                contable << graderrornorm;
                contable << residualnorm;
                contable << Float( end - start );
                contable << nl;
                
                contable.lg();

                {
                    fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                    VTKWriter vtk( M, fs, getbasename(__FILE__) );
                    
                    if( r == 1 ) { 
                        vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                    } 
                    
                    {
                        const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );
                        const auto printable_sol = interpol_matrix * incmatrix * sol; 
                        vtk.writeCellScalarData( printable_sol, "iterativesolution_scalar_data_cellwise" , 1.0 );
                    }
                    
                    {
                        auto computed_grad = aug_diffmatrix * elevmatrix * incmatrix * sol;
                        const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r-1 );
                        const auto printable_grad = interpol_matrix * computed_grad; 
                        vtk.writeCellVectorData_barycentricgradients( printable_grad, "gradient_interpolation" , 1.0 );
                    }
                    
                    fs.close();
                }

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        
            } 
        
        }
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
