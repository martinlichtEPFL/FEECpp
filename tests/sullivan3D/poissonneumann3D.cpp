

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
// #include "../../fem/global.lagrangeincl.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: 3D Poisson Neumann Problem" << nl;
        
        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial3D M = StandardCube3D();
            
            M.check();
            
            LOG << "Prepare scalar fields for testing..." << nl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 3 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            
            const Float xfeq = 1.;
            const Float yfeq = 1.;
            const Float zfeq = 1.;
            

            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    // return FloatVector({ 1. });
                    return FloatVector({ 
                           std::cos( xfeq * Constants::pi * vec[0] )
                         * std::cos( yfeq * Constants::pi * vec[1] )
                         * std::cos( zfeq * Constants::pi * vec[2] ) 
                         });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    return FloatVector( { 
                            -xfeq * Constants::pi * std::sin( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) * std::cos( zfeq * Constants::pi * vec[2] ),
                            -yfeq * Constants::pi * std::cos( xfeq * Constants::pi * vec[0] ) * std::sin( yfeq * Constants::pi * vec[1] ) * std::cos( zfeq * Constants::pi * vec[2] ), 
                            -zfeq * Constants::pi * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) * std::sin( zfeq * Constants::pi * vec[2] ), 
                        });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    return FloatVector({ 
                        xfeq*xfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) * std::cos( yfeq * Constants::pi * vec[2] )
                        +
                        yfeq*yfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) * std::cos( yfeq * Constants::pi * vec[2] )
                        +
                        zfeq*zfeq * Constants::pisquare * std::cos( xfeq * Constants::pi * vec[0] ) * std::cos( yfeq * Constants::pi * vec[1] ) * std::cos( yfeq * Constants::pi * vec[2] )
                     });
                };
            

            

            

            const int min_l = 0; 
            const int max_l = 4;
            
            const int min_r = 1;
            const int max_r = 1;
            
            ConvergenceTable contable("Mass error");
            
            contable << "u_error" << "du_error" << "residual" << "time" << nl;
            

            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << nl;
                LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ ) 
                {
                    
                    LOG << "...assemble scalar mass matrices" << nl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );

                    LOG << "...assemble vector mass matrix" << nl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );
                    
                    LOG << "...assemble differential matrix and transpose" << nl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    LOG << "...assemble inclusion matrix and transpose" << nl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    LOG << "...assemble stiffness matrix" << nl;
            
                    // ProductOperator 
//                     auto stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    auto op1 = incmatrix_t * diffmatrix_t;
                    auto op2 = op1 * vector_massmatrix;
                    auto op3 = op2 * diffmatrix;
                    auto stiffness = op3 * incmatrix;
                    auto& stiffness_csr = stiffness;

//                     auto opr  = diffmatrix & incmatrix;
//                     auto opl  = opr.getTranspose(); 
//                     auto stiffness = opl & ( vector_massmatrix & opr );
                    
//                     stiffness.sortentries();
//                     auto stiffness_csr = MatrixCSR( stiffness );
                    
                    auto stiffness_invprecon = DiagonalOperator( stiffness.getdimin(), 1. );
//                     auto stiffness_invprecon = InverseDiagonalPreconditioner( stiffness );
                    LOG << "Average value of diagonal preconditioner: " << stiffness_invprecon.getdiagonal().average() << nl;

                    {

                        const auto& function_sol  = experiment_sol;
                        const auto& function_grad = experiment_grad;
                        const auto& function_rhs  = experiment_rhs;
                        
                        LOG << "...interpolate explicit solution and rhs" << nl;
            
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        
                        FloatVector interpol_one  = Interpolation( M, M.getinnerdimension(), 0, r, constant_one );
                        
                        LOG << "...measure kernel component: ";
            
                        Float average_sol = interpol_one * ( scalar_massmatrix * interpol_sol );
                        Float average_rhs = interpol_one * ( scalar_massmatrix * interpol_rhs );
                        
                        LOG << average_sol << space << average_rhs << nl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( incmatrix.getdimin(), 0. );
                        
                        LOG << "...iterative solver" << nl;
                        
                        timestamp start = timestampnow();
                        
                        {
                            sol.zero();
                            ConjugateResidualMethod Solver( stiffness_csr );
                            Solver.solve( sol, rhs );
                        }

                        timestamp end = timestampnow();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                        
                        LOG << "...compute error and residual" << nl;
                        
                        auto computed_sol  = incmatrix * sol;
                        auto computed_grad = diffmatrix * incmatrix * sol;

                        auto polluted_error = interpol_sol  - incmatrix * sol;
                        Float average_error = ( interpol_one * ( scalar_massmatrix * polluted_error ) ) / ( interpol_one * ( scalar_massmatrix * interpol_one ) );
                        auto          error = polluted_error - average_error * interpol_one;
                        auto      graderror = interpol_grad - diffmatrix * incmatrix * sol;
                        
                        Float errornorm       = std::sqrt( error * ( scalar_massmatrix * error ) );
                        Float graderrornorm   = std::sqrt( graderror * ( vector_massmatrix * graderror ) );
                        Float residualnorm    = ( rhs - stiffness * sol ).norm();
                        
                        LOG << "error:     " << errornorm    << nl;
                        LOG << "graderror: " << graderrornorm << nl;
                        LOG << "residual:  " << residualnorm << nl;
                        
                        contable << errornorm;
                        contable << graderrornorm;
                        contable << residualnorm;
                        contable << Float( end - start );
                        contable << nl;
                        
                        contable.lg();


                        {
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                            VTKWriter vtk( M, fs, getbasename(__FILE__) );
                            
                            {
                                const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );
                                const auto printable_sol = interpol_matrix * incmatrix * sol; 
                                vtk.writeCellScalarData( printable_sol, "iterativesolution_scalar_data_cellwise" , 1.0 );
                            }
                            
                            {
                                const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r-1 );
                                const auto printable_grad = interpol_matrix * computed_grad; 
                                vtk.writeCellVectorData_barycentricgradients( printable_grad, "gradient_interpolation" , 1.0 );
                            }
                            
                            fs.close();
                        }


                    }
                    
                }

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
                
            } 
        
        }
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
