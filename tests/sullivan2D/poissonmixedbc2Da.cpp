

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: 2D Poisson problem, Dirichlet BC on one side " << nl;
        
        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.set_flag( 1, 0, SimplexFlag::SimplexFlagDirichlet );
            M.set_flag( 0, M.get_subsimplex( 1, 0, 0, 0 ), SimplexFlag::SimplexFlagDirichlet );
            M.set_flag( 0, M.get_subsimplex( 1, 0, 0, 1 ), SimplexFlag::SimplexFlagDirichlet );

            
            LOG << "Prepare scalar fields for testing..." << nl;
            

            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            
            
            
            


            
            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    Float x = vec[0]; Float y = vec[1];
                    return FloatVector({ square((square(x)-1)*(square(y)-1)) });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    Float x = vec[0]; Float y = vec[1];
                    return FloatVector( { 
                        4 * x * ( x*x - 1 ) * square( y*y - 1 ),
                        4 * y * square( x*x - 1 ) * ( y*y - 1 )  
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    Float x  =  vec[0]; Float y  =  vec[1];
                    Float x2 =     x*x; Float y2 =     y*y;
                    Float x4 = x*x*x*x; Float y4 = y*y*y*y;
                    return FloatVector({ 
                        - 4 * ( x4 * ( 3 * y2 - 1 ) + x2 * ( 3 * y4 - 12 * y2 + 5 ) - y4 + 5 * y2 - 2)  
                     });
                };
            

            

            

            const int min_l = 0; 
            const int max_l = 3;
            
            const int min_r = 3;
            const int max_r = 3;
            
            ConvergenceTable contable("Mass error");
            
            contable << "u_error" << "du_error" << "residual" << "time" << nl;
            

            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << nl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                if( l != 0 )
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
//                     auto op1 = incmatrix_t * diffmatrix_t;
//                     auto op2 = op1 * vector_massmatrix;
//                     auto op3 = op2 * diffmatrix;
//                     auto stiffness = op3 * incmatrix;
//                     auto& stiffness_csr = stiffness;

                    auto opr = diffmatrix & incmatrix;
                    auto opl = opr.getTranspose(); 
                    auto stiffness = opl & ( vector_massmatrix & opr );
                    
                    stiffness.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness );
                    
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
                        
                        LOG << "...compute norms of solution and right-hand side" << nl;
            
                        Float sol_norm = interpol_sol * ( scalar_massmatrix * interpol_sol );
                        Float rhs_norm = interpol_rhs * ( scalar_massmatrix * interpol_rhs );
                        
                        LOG << "solution norm: " << sol_norm << nl;
                        LOG << "rhs norm:      " << rhs_norm << nl;

                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( incmatrix.getdimin(), 0. );
                        
                        LOG << "...iterative solver" << nl;
                        

                        timestamp start = timestampnow();

                        {
                            LOG << "CGM - Classic" << nl;
                        
                            sol.zero();
                            
                            FloatVector residual( rhs );
                            
                            const FloatVector diagonal = stiffness_csr.getDiagonal();
                            assert( diagonal.isfinite() );
                            assert( diagonal.isnonnegative() );
                            
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

                        timestamp end = timestampnow();
                        LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                        LOG << "...compute error and residual" << nl;
            
                        
                        auto computed_sol  = incmatrix * sol;
                        auto computed_grad = diffmatrix * incmatrix * sol;
                        
                        auto errornorm_aux = interpol_sol  - computed_sol;
                        auto graderror_aux = interpol_grad - computed_grad;
                        
                        Float errornorm     = sqrt( errornorm_aux * ( scalar_massmatrix * errornorm_aux ) );
                        Float graderrornorm = sqrt( graderror_aux * ( vector_massmatrix * graderror_aux ) );
                        Float residualnorm  = ( rhs - stiffness * sol ).norm();
                        
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
                            
                            if( r == 1 ) { 
                                vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                            } 
                            
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
