

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../utility/math.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsolver.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"
#include "../../vtk/vtkwriter.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test: 2D curl-curl problem" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial2D M = StandardSquare2D_simple();
    
    M.check();
    
    M.automatic_dirichlet_flags();
    M.check_dirichlet_flags();

    
    LOG << "Prepare scalar fields for testing..." << nl;
    

    std::function<FloatVector(const FloatVector&)> constant_one
        = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ 1. });
            };
    
    
    // u dx + v dy -> ( v_x - u_y ) dxdy -> ( - u_yy + v_xy , - v_xx + u_xy )
    
    // const Float A = Constants::twopi;

            
    std::function<FloatVector(const FloatVector&)> experiment_sol = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 2 );
            // return FloatVector({ 1. });
            return FloatVector({ 
                 //blob( vec[0] ) * blob_dev( vec[1] )
                 square( vec[0]*vec[0] - 1. ) * 4. * vec[1] * ( vec[1]*vec[1] - 1. )
                ,
                //-blob_dev( vec[0] ) * blob( vec[1] )
                -4. * vec[0] * ( vec[0]*vec[0] - 1. ) * square( vec[1]*vec[1] - 1. )
                });
        };
    
    std::function<FloatVector(const FloatVector&)> experiment_curl = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 2 );
            return FloatVector({ 17. });
        };
    
    std::function<FloatVector(const FloatVector&)> experiment_rhs = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 2 );
            return FloatVector({ 
                blob_dev(vec[0])*blob(vec[1]) 
                + 
                // ( - blob( vec[0] ) * blob_devdevdev( vec[1] )     - blob_devdev( vec[0] ) * blob_dev( vec[1] ) )
                ( - square( vec[0]*vec[0] - 1. ) * 24. * vec[1]      - ( 12. * vec[0]*vec[0] - 4. ) * ( 4. * vec[1] * ( vec[1]*vec[1] - 1. ) ) )
                ,
                blob(vec[0])*blob_dev(vec[1])
                + 
                // ( + blob_devdevdev( vec[0] ) * blob_dev( vec[1] ) + blob_dev( vec[0] ) * blob_devdev( vec[1] ) )
                (   24. * vec[0] * square( vec[1]*vec[1] - 1. )      + ( 4. * vec[0] * ( vec[0]*vec[0] - 1. ) ) * ( 12. * vec[1]*vec[1] - 4. ) )
                });
        };
    
    std::function<FloatVector(const FloatVector&)> experiment_aux = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 2 );
            // return FloatVector({ 1. });
            return FloatVector( 1, blob(vec[0])*blob(vec[1]) );
        };
    
    
    
    
    

    

    const int min_l = 1; 
    const int max_l = 4;
    
    const int min_r = 1;
    const int max_r = 1;
    
    
    ConvergenceTable contable("Mass error and numerical residuals");
    
    contable << "u_error" << "du_error" << "sigma_error" << "u_res" << "sigma_res" << "time" << nl;

    
    assert( 0 <= min_l and min_l <= max_l );
    assert( 0 <= min_r and min_r <= max_r );
        
    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        for( int r = min_r; r <= max_r; r++ ) 
        {
            
            LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
            LOG << "... assemble matrices" << nl;
    
            SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

            LOG << "... assemble inclusion matrices" << nl;
    
            SparseMatrix scalar_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

            LOG << "... assemble algebraic matrices" << nl;
    
            SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r+1 );
            SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

            SparseMatrix vector_incmatrix   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r   );
            SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

            SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
            SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();

            LOG << "... compose system matrices" << nl;
    
            auto mat_A  = vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix;
            mat_A.sortandcompressentries();
                
            auto mat_Bt = vector_incmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix; // upper right
            mat_Bt.sortandcompressentries();
            
            auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
            mat_B.sortandcompressentries();
            
            LOG << "... compose CSR system matrices" << nl;
    
            auto A  = MatrixCSR( mat_A  );
            auto Bt = MatrixCSR( mat_Bt );
            auto B  = MatrixCSR( mat_B  );
            
            auto C  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
            
            // TODO: develop preconditioners 
            // auto PA = IdentityMatrix( A.getdimin() );
            // auto PC = IdentityMatrix( C.getdimin() );
            auto PA = MatrixCSR( vector_incmatrix_t & vector_massmatrix & vector_incmatrix )
                              + MatrixCSR( vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix );
            auto PC = MatrixCSR( scalar_incmatrix_t & scalar_diffmatrix_t & vector_massmatrix & scalar_diffmatrix & scalar_incmatrix );
                
            LOG << "share zero PA = " << PA.getnumberofzeroentries() << "/" <<  PA.getnumberofentries() << nl;
            LOG << "share zero PC = " << PC.getnumberofzeroentries() << "/" <<  PC.getnumberofentries() << nl;
                        
                        
            const auto& function_sol = experiment_sol;
            const auto& function_rhs = experiment_rhs;
            const auto& function_aux = experiment_aux;
            const auto& function_curl = experiment_curl;
            
            LOG << "...interpolate explicit solution and rhs" << nl;
            
            FloatVector interpol_sol  = Interpolation( M, M.getinnerdimension(), 1, r,    function_sol );
            FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 1, r,    function_rhs );
            FloatVector interpol_aux  = Interpolation( M, M.getinnerdimension(), 0, r+1,  function_aux );
            FloatVector interpol_curl = Interpolation( M, M.getinnerdimension(), 2, r-1, function_curl );
            
            FloatVector rhs_sol = vector_incmatrix_t * vector_massmatrix * interpol_rhs;
            FloatVector rhs_aux = scalar_incmatrix_t * scalar_diffmatrix_t * vector_massmatrix * interpol_sol ;// FloatVector( B.getdimout(), 0. );

            FloatVector sol( A.getdimout(), 0. );
            FloatVector aux( B.getdimout(), 0. );

            // compute the solution ....



            timestamp start = timestampnow();

            //TODO: set up operator preconditioner
            // { 
            //     auto X = Block2x2Operator( A.getdimout() + B.getdimout(), A.getdimin() + Bt.getdimin(), A, Bt, B, C );    
            //     //
            //     FloatVector sol_full( A.getdimin()  + Bt.getdimin(),  0. );
            //     sol_full.setslice(             0, sol );
            //     sol_full.setslice( A.getdimout(), aux );
            //     //
            //     FloatVector rhs_full( A.getdimout() +  B.getdimout(), 0. );
            //     rhs_full.setslice(             0, rhs_sol );
            //     rhs_full.setslice( A.getdimout(), rhs_aux );
            //     //
            //     HerzogSoodhalterMethod Solver( X );
            //     Solver.solve( sol_full, rhs_full );
            //     //
            //     sol = sol_full.getslice(             0, A.getdimout() );
            //     aux = sol_full.getslice( A.getdimout(), B.getdimout() );
            // }

            {
                const auto PAinv = inv(PA,desired_precision,-1);
                const auto PCinv = inv(PC,desired_precision,-1);
                BlockHerzogSoodhalterMethod( 
                    sol, 
                    aux, 
                    rhs_sol, 
                    rhs_aux, 
                    A, Bt, B, C, 
                    desired_precision,
                    1,
                    PAinv, PCinv
                );
            }

            timestamp end = timestampnow();

            // ... computed the solution

            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
            
            LOG << "...compute error and residual" << nl;

            auto errornorm_aux_sol  = interpol_sol  - vector_incmatrix * sol;
            auto errornorm_aux_curl = interpol_curl - vector_diffmatrix * vector_incmatrix * sol;
            auto errornorm_aux_aux  = interpol_aux  - scalar_incmatrix * aux;

            Float errornorm_sol  = sqrt( errornorm_aux_sol  * ( vector_massmatrix * errornorm_aux_sol ) );
            Float errornorm_curl = sqrt( errornorm_aux_curl * ( volume_massmatrix * errornorm_aux_curl ) );
            Float errornorm_aux  = sqrt( errornorm_aux_aux  * ( scalar_massmatrix * errornorm_aux_aux ) );
            
            Float residual_sol  = ( rhs_sol - A * sol - Bt * aux ).norm();
            Float residual_aux  = ( rhs_aux - B * sol -  C * aux ).norm();

            LOG << "error:        " << errornorm_sol << nl;
            LOG << "aux error:    " << errornorm_aux << nl;
            LOG << "residual:     " << residual_sol << nl;
            LOG << "aux residual: " << residual_aux << nl;

            contable << errornorm_sol;
            contable << errornorm_curl;
            contable << errornorm_aux;
            contable << residual_sol;
            contable << residual_aux;
            contable << Float( end - start );
            contable << nl;

            contable.lg();
            


            {
                fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                VTKWriter vtk( M, fs, getbasename(__FILE__) );

                {
                    const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 0, 0, r+1 );
                    const auto printable_aux = interpol_matrix * scalar_incmatrix * aux; 
                    vtk.writeCellScalarData( printable_aux, "computed_aux" , 1.0 );
                }
                
                {
                    const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );
                    const auto printable_sol = interpol_matrix * vector_incmatrix * sol; 
                    vtk.writeCellVectorData_barycentricgradients( printable_sol, "computed_solution" , 1.0 );
                }
                
                {
                     const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 2, 0, r-1 );
                     const auto printable_curl = interpol_matrix * vector_diffmatrix * vector_incmatrix * sol; 
                     Assert( printable_curl.getdimension() == (M.getinnerdimension()+1) * M.count_simplices(M.getinnerdimension()), 
                                    printable_curl.getdimension(), M.count_simplices(M.getinnerdimension()) );
                     vtk.writeCellScalarData_barycentricvolumes( printable_curl, "computed_curl" , 1.0 );
                }
                            
                assert( function_aux( FloatVector{0.0,0.0 }).getdimension() == 1 );
                vtk.writeCellScalarData( function_aux,  "function_aux" , 1.0 );
                vtk.writeCellVectorData( function_sol,  "function_sol" , 1.0 );
                vtk.writeCellScalarData( function_curl, "function_curl" , 1.0 );

                vtk.writeCellVectorData( function_rhs,  "function_rhs" , 1.0 );
            
                fs.close();
            }
            
        }

        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

    } 
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
