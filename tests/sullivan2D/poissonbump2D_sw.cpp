

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
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/global.elevation.hpp"
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
            const int max_l = 4;
            
            const int min_r = 1;
            const int max_r = 1;

            
            const int sullivan = 0;
            const int whitney  = 1;
            ConvergenceTable contable[2] = { ConvergenceTable("Mass error"), ConvergenceTable("Mass error") };
            
            contable[sullivan] << "u_error" << "du_error" << "residual" << "time" << nl;
            contable[whitney]  << "u_error" << "du_error" << "residual" << "time" << nl;
            

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
            
                    SparseMatrix incmatrix_s = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    
                    SparseMatrix incmatrix_s_t = incmatrix_s.getTranspose();

                    SparseMatrix incmatrix_w = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r );
                    // FEECBrokenElevationMatrix( M, M.getinnerdimension(), 0, r, 1 ) & 
                    
                    
                    SparseMatrix incmatrix_w_t = incmatrix_w.getTranspose();

                    LOG << "...assemble stiffness matrix" << nl;
            
                    auto opr_s  = diffmatrix & incmatrix_s;
                    auto opl_s  = opr_s.getTranspose(); 
                    auto stiffness_s = opl_s & ( vector_massmatrix & opr_s );
                    
                    stiffness_s.sortentries();
                    auto stiffness_csr_s = MatrixCSR( stiffness_s );
                    
                    auto opr_w  = diffmatrix & incmatrix_w;
                    auto opl_w  = opr_w.getTranspose(); 
                    auto stiffness_w = opl_w & ( vector_massmatrix & opr_w );
                    
                    stiffness_w.sortentries();
                    auto stiffness_csr_w = MatrixCSR( stiffness_w );
                    
                    
                    auto stiffness_invprecon_s = DiagonalOperator( stiffness_s.getdimin(), 1. );
                    auto stiffness_invprecon_w = DiagonalOperator( stiffness_w.getdimin(), 1. );
                    LOG << "Average value of diagonal preconditioner (Sullivan): " << stiffness_invprecon_s.getdiagonal().average() << nl;
                    LOG << "Average value of diagonal preconditioner (Whitney) : " << stiffness_invprecon_w.getdiagonal().average() << nl;

                    {

                        const auto& function_sol  = experiment_sol;
                        const auto& function_grad = experiment_grad;
                        const auto& function_rhs  = experiment_rhs;
                        
                        LOG << "...interpolate explicit solution and rhs" << nl;
            
                        FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, r,   function_sol  );
                        FloatVector interpol_grad = Interpolation( M, M.getinnerdimension(), 1, r-1, function_grad );
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        
                        FloatVector rhs_s = incmatrix_s_t * ( scalar_massmatrix * interpol_rhs );
                        FloatVector sol_s( incmatrix_s.getdimin(), 0. );
                        
                        FloatVector rhs_w = incmatrix_w_t * ( scalar_massmatrix * interpol_rhs );
                        FloatVector sol_w( incmatrix_w.getdimin(), 0. );
                        

                        LOG << "...iterative solver (Sullivan)" << nl;
                        timestamp start_s = timestampnow();

                        {
                            sol_s.zero();
                            ConjugateGradientMethod Solver( stiffness_csr_s );
                            Solver.solve( sol_s, rhs_s );
                        }

                        timestamp end_s = timestampnow();
                        LOG << "\t\t\t Time (Sullivan): " << timestamp2measurement( end_s - start_s ) << nl;

                        LOG << "...iterative solver (Whitney)" << nl;
                        timestamp start_w = timestampnow();

                        {
                            sol_w.zero();
                            ConjugateGradientMethod Solver( stiffness_csr_w );
                            Solver.solve( sol_w, rhs_w );
                        }

                        timestamp end_w = timestampnow();
                        LOG << "\t\t\t Time (Whitney): " << timestamp2measurement( end_w - start_w ) << nl;



                        LOG << "...compute error and residual (Sullivan)" << nl;
                        
                        auto computed_sol_s  = incmatrix_s * sol_s;
                        auto computed_grad_s = diffmatrix * incmatrix_s * sol_s;
                        
                        auto errornorm_aux_s = interpol_sol  - computed_sol_s;
                        auto graderror_aux_s = interpol_grad - computed_grad_s;
                        
                        Float errornorm_s     = sqrt( errornorm_aux_s * ( scalar_massmatrix * errornorm_aux_s ) );
                        Float graderrornorm_s = sqrt( graderror_aux_s * ( vector_massmatrix * graderror_aux_s ) );
                        Float residualnorm_s  = ( rhs_s - stiffness_s * sol_s ).norm();
                        
                        LOG << "error:     " << errornorm_s     << nl;
                        LOG << "graderror: " << graderrornorm_s << nl;
                        LOG << "residual:  " << residualnorm_s  << nl;
                        LOG << "time:      " << Float( end_s - start_s ) << nl;
                        
                        contable[sullivan] << errornorm_s;
                        contable[sullivan] << graderrornorm_s;
                        contable[sullivan] << residualnorm_s;
                        contable[sullivan] << Float( end_s - start_s );
                        contable[sullivan] << nl;
                        
                        contable[sullivan].lg();



                        LOG << "...compute error and residual (Whitney)" << nl;
                        
                        auto computed_sol_w  = incmatrix_w * sol_w;
                        auto computed_grad_w = diffmatrix * incmatrix_w * sol_w;
                        
                        auto errornorm_aux_w = interpol_sol  - computed_sol_w;
                        auto graderror_aux_w = interpol_grad - computed_grad_w;
                        
                        Float errornorm_w     = sqrt( errornorm_aux_w * ( scalar_massmatrix * errornorm_aux_w ) );
                        Float graderrornorm_w = sqrt( graderror_aux_w * ( vector_massmatrix * graderror_aux_w ) );
                        Float residualnorm_w  = ( rhs_w - stiffness_w * sol_w ).norm();
                        
                        LOG << "error:     " << errornorm_w     << nl;
                        LOG << "graderror: " << graderrornorm_w << nl;
                        LOG << "residual:  " << residualnorm_w  << nl;
                        LOG << "time:      " << Float( end_w - start_w ) << nl;
                        
                        contable[whitney] << errornorm_w;
                        contable[whitney] << graderrornorm_w;
                        contable[whitney] << residualnorm_w;
                        contable[whitney] << Float( end_w - start_w );
                        contable[whitney] << nl;
                        
                        contable[whitney].lg();


                        {
                            
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                            VTKWriter vtk( M, fs, getbasename(__FILE__) );
                            
                            
                            if( r == 1) {
                                vtk.writeVertexScalarData( sol_s, "iterativesolution_scalar_data_s" , 1.0 );
                                vtk.writeVertexScalarData( sol_w, "iterativesolution_scalar_data_w" , 1.0 );
                            }
                            
                            {
                                const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r );
                                const auto printable_sol = interpol_matrix * incmatrix_s * sol_s; 
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
                                const auto printable_grad = interpol_matrix * computed_grad_s; 
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
