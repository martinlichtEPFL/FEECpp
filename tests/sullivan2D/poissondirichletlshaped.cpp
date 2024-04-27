

/**/

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../dense/factorization.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../sparse/rainbow.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../fem/lagrangematrices.hpp"
#include "../../fem/utilities.hpp"



FloatVector IncreaseResolution( const MeshSimplicial2D& mesh, const FloatVector& lores );

FloatVector IncreaseResolution( const MeshSimplicial2D& mesh, const FloatVector& lores )
{
    int counter_vertices = mesh.count_vertices();
    int counter_edges    = mesh.count_edges();

    assert( lores.getdimension() == counter_vertices );

    FloatVector hires( counter_vertices + counter_edges );

    for( int v = 0; v < counter_vertices; v++ ) hires[v] = lores[v];

    for( int e = 0; e < counter_edges; e++ ) {
        int v0 = mesh.get_edge_vertex(e,0);
        int v1 = mesh.get_edge_vertex(e,1);
        Float value = ( lores[v0] + lores[v1] ) / 2.;
        hires[counter_vertices + e] = value;
    }

    return hires;
}







using namespace std;

int main( int argc, char *argv[] )
{
        
    LOG << "Unit Test: 2D Dirichlet problem on L-shaped domain" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial2D M = LShapedDomain2D();

    M.automatic_dirichlet_flags();
    
    M.check();
    
    M.check_dirichlet_flags();
    
    LOG << "Prepare scalar fields for testing..." << nl;
    

    std::function<FloatVector(const FloatVector&)> constant_one
        = [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 2 );
                return FloatVector({ 1. });
            };
    
    
    
    


    // The solution of Laplacian problems over L-shaped domainswith a singular function boundary integral method
    // https://onlinelibrary.wiley.com/doi/pdf/10.1002/cnm.489?casa_token=KTbdSboKSK8AAAAA:ISbMXTrwR6i-CocYB6hgQdxdbGgjQxo1QMxRA-L97XFrW_BuEiyUxnXZSVM_SF3DTLmHyGe0ZNdLZtR3
            
    std::function<FloatVector(const FloatVector&)> experiment_rhs = 
        [=](const FloatVector& vec) -> FloatVector{
            assert( vec.getdimension() == 2 );
            return FloatVector({
                1.0
                });
        };
    

    

        

    const int min_l = 1; 
    const int max_l = 6;

    std::vector<MeshSimplicial2D>    meshes;
    std::vector<FloatVector>         solutions;
    std::vector<ConvergenceTable>    contables;
    
    ConvergenceTable contable_pastone("Mass errror, one before");
    ConvergenceTable contable_pasttwo("Mass errror, two before");   

    contable_pastone << "level" << "u_error" << "du_error" << nl;
    contable_pasttwo << "level" << "u_error" << "du_error" << nl;


    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        {
            
            const int r = 1;
            const int w = r;

            LOGPRINTF("Polynomial degrees: r=%d w=%d \n", r, w );
            
            LOG << "...assemble scalar mass matrices" << nl;
    
            auto mass = MatrixCSR( LagrangeMassMatrix( M, r ) );

            LOG << "...assemble vector mass matrix" << nl;
    
            auto stiffness = MatrixCSR( LagrangeStiffnessMatrix( M, r ) );
            
            mass.compressentries();
            stiffness.compressentries();

            LOG << "Size of mass matrix:      " << mass.memorysize() << nl;
            LOG << "Size of stiffness matrix: " << stiffness.memorysize() << nl;

            assert( mass.isfinite() );
            assert( stiffness.isfinite() );



            {

                LOG << "...interpolate explicit solution and rhs" << nl;
    
                FloatVector      sol( M.count_vertices(), 0. );
                
                FloatVector residual( M.count_vertices(), 0. );

                FloatVector     rhs( M.count_vertices(), 0. );
                
                #if defined(_OPENMP)
                #pragma omp parallel for
                #endif
                for( int s = 0; s < M.count_triangles(); s++ )
                {
                    auto measure = M.getMeasure(2,s);

                    auto midpoint = M.get_midpoint(2,s);

                    auto f = experiment_rhs(midpoint);

                    auto contribution = f[0] * measure / 3.;

                    for( int vi = 0; vi <= 2; vi++ )
                    {
                        int v = M.get_subsimplex( 2, 0, s, vi );
                        
                        if( M.get_flag(0,v) == SimplexFlag::SimplexFlagDirichlet ) continue;

                        #if defined(_OPENMP)
                        #pragma omp atomic
                        #endif
                        rhs[v] += contribution;

                    }

                }

                assert( rhs.isfinite() );

                
                LOG << "...iterative solver" << nl;
                
                auto& stiffness_csr = stiffness;

                sol.zero();

                auto diagonal = stiffness.getDiagonal();
            
                Rainbow rainbow( stiffness );

                timestamp solver_time;

                timestamp start = timestampnow();
                ConjugateGradientSolverCSR_Rainbow( 
                    sol.getdimension(), 
                    sol.raw(), 
                    rhs.raw(), 
                    stiffness_csr.getA(), stiffness_csr.getC(), stiffness_csr.getV(),
                    residual.raw(),
                    desired_precision,
                    0,
                    diagonal.raw(),
                    1.2 // empircally chosen at 1.2, barely any influence
                    , rainbow.num_colors, rainbow.F.data(), rainbow.B.data(), rainbow.R.data()
                );
                timestamp end = timestampnow();
                solver_time = end - start;
                LOG << "\t\t\t Time: " << timestamp2measurement( solver_time ) << nl;

                residual = rhs - stiffness * sol;

                Float residualnorm = residual.norm();

                LOG << "residual:  " << residualnorm << nl;
                LOG << "time:      " << Float( solver_time ) << nl;
                
                
                LOG << "...update saved old solutions:" << nl;
                if( l > min_l )
                {
                    std::vector<FloatVector> new_solutions;

                    const auto& old_M = meshes[ l - min_l - 1 ];

                    for( const auto& sol : solutions )
                        new_solutions.push_back( IncreaseResolution( old_M, sol ) );

                    assert( solutions.size() == new_solutions.size() );
                    
                    solutions.clear();

                    assert( solutions.size() == 0 );

                    solutions = std::move( new_solutions );
                }

                solutions.push_back( sol );

                Assert( solutions.size() == l-min_l+1, solutions.size(), l-min_l+1 );



                LOG << "...compute errors against previous solutions:" << nl;
                
                ConvergenceTable contable( printf_into_string("Mass error (l=%d)", l ) );
                
                contable << "level" << "u_error" << "du_error" << nl;

                for( int i = 0; i < solutions.size(); i++ )
                {

                    const auto& old_sol = solutions[i];

                    FloatVector error   = sol - old_sol;

                    assert( error.isfinite() );
                    assert( ( mass * error ).isfinite() );



                    LOG <<  error * ( mass * error ) << nl;
                    LOG <<  error * ( stiffness * error ) << nl;
                    
                    Float errornorm     = std::sqrt( error * ( mass * error ) );
                    Float graderrornorm = std::sqrt( error * ( stiffness * error ) );
                    
                    LOG << "error:     " << errornorm    << nl;
                    LOG << "graderror: " << graderrornorm << nl;
                            
                    contable << Float(i+min_l);
                    contable << errornorm;
                    contable << graderrornorm;
                    contable << nl;

                    if( i+min_l == l-1 )
                        contable_pastone << Float(l-1) << errornorm << graderrornorm << nl;
                    
                    if( i+min_l == l-2 )
                        contable_pasttwo << Float(l-2) << errornorm << graderrornorm << nl;
                    
                }

                contable.lg();


                // save the mesh, solutions, and convergence tables

                meshes.push_back(M);
                // solutions.push_back(sol);
                contables.push_back(contable);

                
                {
                    fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                    VTKWriter vtk( M, fs, getbasename(__FILE__) );
                    
                    assert( r == 1 );

                    vtk.writeVertexScalarData( sol, "iterativesolution_scalar_data" , 1.0 );
                    
                    fs.close();
                }


            }
            
        }

        for( const auto& contable: contables )
            contable.lg();

        contable_pastone.lg();
        contable_pasttwo.lg();
        
        if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }
        
        

    } 
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
