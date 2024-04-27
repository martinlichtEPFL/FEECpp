

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/math.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: 2D Dirichlet problem with Bump function" << nl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
           
            M.check_dirichlet_flags();

            M.getcoordinates().scale(1.1);
            
            LOG << "Prepare scalar fields for testing..." << nl;
            


            std::function<FloatVector(const FloatVector&)> experiment_sol = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector({ 
                        bumpfunction(vec[0]) * bumpfunction(vec[1])
                    });
                };
            
            std::function<FloatVector(const FloatVector&)> experiment_grad = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    // return FloatVector({ 1. });
                    return FloatVector( { 
                            bumpfunction_dev(vec[0]) *     bumpfunction(vec[1]),
                            bumpfunction(vec[0])     * bumpfunction_dev(vec[1]), 
                    });
                };
            

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        -
                        bumpfunction_devdev(vec[0]) *        bumpfunction(vec[1])
                        -
                        bumpfunction(vec[0])        * bumpfunction_devdev(vec[1])
                    });
                };
            
            
            

            

            const int min_l = 0; 
            const int max_l = 5;
            
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
            
                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
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
                        
                        FloatVector rhs = incmatrix_t * ( scalar_massmatrix * interpol_rhs );

                        FloatVector sol( incmatrix.getdimin(), 0. );
                        
                        LOG << "...iterative solver" << nl;
                        
                        timestamp start = timestampnow();

                        {
                            sol.zero();
                            ConjugateGradientMethod Solver( stiffness_csr );
                            Solver.solve( sol, rhs );
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
                            
                            
                            if( r == 1) {
                                vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                            }
                            
                            {
                                const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );
                                const auto printable_sol = interpol_matrix * incmatrix * sol; 
                                vtk.writeCellScalarData( printable_sol, "iterativesolution_scalar_data_cellwise" , 1.0 );
                            }
                            
                            {
                                auto rhs = FloatVector( M.count_simplices(0) );
                                
                                for( int c = 0; c < M.count_simplices(0); c++ ) { 
                                    auto x = M.getcoordinates().getdata(c,0);
                                    auto y = M.getcoordinates().getdata(c,1);
                                    auto value = experiment_rhs( FloatVector({ x, y }) )[0];
                                    rhs[c] = value;
                                }
                                    
                                vtk.writeVertexScalarData( rhs, "reference_scalar_data" , 1.0 );
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