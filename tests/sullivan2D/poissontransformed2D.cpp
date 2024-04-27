

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/global.coefficientmassmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: 2D Poisson transformed problem" << nl;
        
        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D M = RhombicAnnulus2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
           
            M.check_dirichlet_flags();

            
            LOG << "Prepare scalar fields for testing..." << nl;
            

            auto trafo     = [](const FloatVector& vec) -> FloatVector { 
                assert( vec.getdimension() == 2 );
                return vec.sumnorm() / vec.l2norm() * vec;
            };

            // UNUSED auto trafo_inv = [](const FloatVector& vec) -> FloatVector { 
            //     assert( vec.getdimension() == 2 );
            //     return vec.l2norm() / vec.sumnorm() * vec;
            // };

            // auto trafo_jacobian
            auto jacobian = [](const FloatVector& vec) -> DenseMatrix { 
                assert( vec.getdimension() == 2 );
                
                Float x = vec[0];
                Float y = vec[1];
                Float sx = sign(x);
                Float sy = sign(y);
                
                Float l1  = vec.sumnorm();
                Float l2  = vec.l2norm();
                Float l2c = l2*l2*l2;
                
                return DenseMatrix( 2, 2, {
                    x * sx / l2 - x*x * l1/l2c + l1/l2, x * sy / l2 - y*x * l1/l2c,
                    y * sx / l2 - x*y * l1/l2c,         y * sy / l2 - y*y * l1/l2c + l1/l2 
                });
            };



            
            auto physical_u = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    
                    Float r2 = vec.norm_sq();
                    
                    Float val = 0.25 + 3 * log(r2) / ( 32 * log(2) ) - r2/4.;

                    return FloatVector({ val });
                };
            
            auto physical_gradu = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    Float l2sq = vec.norm_sq();
                    return FloatVector( { 
                            vec[0] * ( 3. / ( log(2) * l2sq ) - 8. ) / 16., 
                            vec[1] * ( 3. / ( log(2) * l2sq ) - 8. ) / 16. // TODO: actual solution gradient
                        });
                };
            
            auto physical_f = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        1.
                     });
                };
            
            auto parametric_u = 
                [=](const FloatVector& vec) -> FloatVector{
                    return physical_u ( trafo( vec ) );
                };
            
            auto parametric_gradu = 
                [=](const FloatVector& vec) -> FloatVector{
                    return Transpose(jacobian(vec)) * physical_gradu ( trafo( vec ) );
                };
            
            auto parametric_f = 
                [=](const FloatVector& vec) -> FloatVector{
                    return physical_f( trafo(vec) );
                };
            
            
            // std::function<DenseMatrix(const FloatVector&)> 
            auto weight_scalar = 
                [=](const FloatVector& vec) -> DenseMatrix{
                    assert( vec.getdimension() == 2 );
                    // return DenseMatrix(1,1, kronecker<int> );
                    return DenseMatrix(1,1,absolute(Determinant(jacobian(vec))));
                };
            
            
            // std::function<DenseMatrix(const FloatVector&)> 
            auto weight_vector = 
                [=](const FloatVector& vec) -> DenseMatrix{
                    assert( vec.getdimension() == 2 );
                    // return DenseMatrix(2,2, kronecker<int> );
                    auto jac = jacobian(vec);

                    auto det = Determinant(jac);

                    // Float x = vec[0];
                    // Float y = vec[1];
                    // Float sx = sign(x);
                    // Float sy = sign(y);
                    // Float other_det = ( ( abs(x) + abs(y) ) * (sx*x+sy*y) / (x*x+y*y) );
                    // LOG << x << tab << y << tab << sx << tab << sy << tab << det << tab << other_det << nl;
                    // assert(
                    //     is_numerically_close( det, ( ( abs(x) + abs(y) ) * (sx*x+sy*y) / (x*x+y*y) ) )
                    // );

                    return absolute(det) * Inverse( Transpose(jac) * jac );
                };
            
            
            
            

            

            const int min_l = 0; 
            const int max_l = 3;

            const int w = 4;
            
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
                
                // if( l != 0 )
                for( int r = min_r; r <= max_r; r++ ) 
                {
                    
                    LOG << "...assemble scalar mass matrices" << nl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 0, r,   w, weight_scalar );

                    LOG << "...assemble vector mass matrix" << nl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenCoefficientMassMatrix( M, M.getinnerdimension(), 1, r-1, w, weight_vector );
                    
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

                        const auto& function_sol  = parametric_u;
                        const auto& function_grad = parametric_gradu;
                        const auto& function_rhs  = parametric_f;
                        
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


                        if( r == 1 ){
                            
                            FloatVector low_interpol_sol  = Interpolation( M, M.getinnerdimension(), 0, 0,   function_sol  );
                            
                            fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                            VTKWriter vtk( M, fs, getbasename(__FILE__) );
                            
                            vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                            vtk.writeCellScalarData( low_interpol_sol, "interpolated_solution" , 1.0 );
                            vtk.writeCellVectorData_barycentricgradients( computed_grad, "gradient_interpolation" , 1.0 );
                            
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
