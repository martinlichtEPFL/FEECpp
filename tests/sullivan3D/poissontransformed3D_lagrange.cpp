

/**/

#include <cmath>

// #include <iostream>

#include <algorithm>
#include <fstream>
#include <vector>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../utility/files.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../sparse/rainbow.hpp"
#include "../../fem/global.coefficientmassmatrix.hpp"
#include "../../fem/lagrangematrices.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"



Float alpha(const Float x);

Float alpha_dev(const Float x);

FloatVector trafo(const FloatVector& vec);

DenseMatrix jacobian(const FloatVector& vec);
        
Float physical_f(const FloatVector& vec) ;

Float parametric_f(const FloatVector& vec);

DenseMatrix physical_A(const FloatVector& vec);

DenseMatrix parametric_A(const FloatVector& vec);

Float weight_scalar(const FloatVector& vec);
        
DenseMatrix weight_vector(const FloatVector& vec);        
        
FloatVector IncreaseResolution( const MeshSimplicial3D& mesh, const FloatVector& lores );





Float alpha(const Float x) { 
    Float maxnorm = absolute(x);
    return std::exp( 1 - 1. / maxnorm );
}

Float alpha_dev(const Float x) { 
    Float maxnorm = absolute(x);
    Float maxnorm_sq = square(maxnorm);            
    return std::exp( 1 - 1. / maxnorm ) / ( maxnorm_sq );
}

const Float B = 2.;

FloatVector trafo(const FloatVector& vec) { 
    assert( vec.getdimension() == 3 );
    
    return vec / B + alpha( vec.maxnorm() ) * ( 1./vec.l2norm() - 1./vec.maxnorm()/B ) * vec;
}

        
DenseMatrix jacobian(const FloatVector& vec) { 
    assert( vec.getdimension() == 3 );
    
    // return DenseMatrix( 3, 3, kronecker<int> );

    Float x = vec[0];
    Float y = vec[1];
    Float z = vec[2];

    Float ax = absolute(x);
    Float ay = absolute(y);
    Float az = absolute(z);

    Float sx = sign(x);
    Float sy = sign(y);
    Float sz = sign(z);

    Float dx = ( ax >= ay and ax >= az ) ? sx : 0.;
    Float dy = ( ay >= ax and ay >= az ) ? sy : 0.;
    Float dz = ( az >= ax and az >= ay ) ? sz : 0.;
    
    Float li  = maximum(ax,maximum(ay,az));
    Float lis = li*li;
    Float l2s = x*x + y*y + z*z;
    Float l2  = sqrt( l2s );
    Float l2c = l2*l2s;

    Float a     = std::exp( 1. - 1. / li );
    Float a_dev = a / lis;
    
    return DenseMatrix( 3, 3, {
        1/B + a * ( 1/l2 - 1/li/B ) + ( a_dev * sx * dx * ( 1/l2 - 1/li/B ) + a * ( -x/l2c - sx*dx/lis/B ) ) * x 
        ,
                                       (                                       a * ( -y/l2c               ) ) * x 
        ,
                                       (                                       a * ( -z/l2c               ) ) * x 
        ,
        //
                                       (                                       a * ( -x/l2c               ) ) * y 
        ,
        1/B + a * ( 1/l2 - 1/li/B ) + ( a_dev * sy * dy * ( 1/l2 - 1/li/B ) + a * ( -y/l2c - sy*dy/lis/B ) ) * y 
        ,
                                       (                                       a * ( -z/l2c               ) ) * y 
        ,
        //
                                       (                                       a * ( -x/l2c               ) ) * z 
        ,
                                       (                                       a * ( -y/l2c               ) ) * z 
        ,
        1/B + a * ( 1/l2 - 1/li/B ) + ( a_dev * sz * dz * ( 1/l2 - 1/li/B ) + a * ( -z/l2c - sz*dz/lis/B ) ) * z 
        ,
    });
}

        
Float physical_f(const FloatVector& vec) 
{
    assert( vec.getdimension() == 3 );
    return 1.;
}

Float parametric_f(const FloatVector& vec) 
{
    return absolute(Determinant(jacobian(vec))) * physical_f( trafo(vec) );
}


DenseMatrix physical_A(const FloatVector& vec) 
{
    assert( vec.getdimension() == 3 );
    Float scalar = vec.l2norm() > 0.5 ? 1. : 2. ;
    return scalar * DenseMatrix(3,3, kronecker<int> );
}

DenseMatrix parametric_A(const FloatVector& vec) 
{
    return physical_A( trafo(vec) );
}





Float weight_scalar(const FloatVector& vec) 
{
    assert( vec.getdimension() == 3 );
    return 1.;
    // return absolute(Determinant(jacobian(vec)));
}
        
DenseMatrix weight_vector(const FloatVector& vec)
{
    assert( vec.getdimension() == 3 );
    
    auto jac = jacobian(vec);
    auto det = Determinant(jac);
    auto invjac = Inverse( jac );
    auto invtjac = Transpose(invjac);
    return absolute(det) * invjac * parametric_A(vec) * Transpose(invjac);

    // auto jac = jacobian(vec);
    // auto det = Determinant(jac);
    // auto invjac = Inverse( jac );
    // auto ret = MatrixTripleMult( parametric_A(vec), invjac );
    // ret *= absolute(det);
    // TransposeSquareInSitu(ret);
    // return ret;
}
        
        


FloatVector IncreaseResolution( const MeshSimplicial3D& mesh, const FloatVector& lores )
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
        
    LOG << "Unit Test for transformed 3D Poisson Problem" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    // MeshSimplicial3D M = FicheraCorner3D();
    MeshSimplicial3D M = CrossedBricks_Five3D();
    
    M.check();

    
    for( int e = 0; e < M.count_edges(); e++ ) {
        Float D = M.getDiameter(1,e);
        if( 1.9 < D and D < 2.1 )
            M.bisect_edge(e); 
    }

    
    for( int s = 0; s < M.count_simplices(2); s++ ) 
    {
        if( M.getsupersimplices(3,2,s).size() > 1 ) continue;
        
        auto midpoint = M.get_midpoint(2,s);

        bool eligible = ( midpoint[1] < 0. );

        if( not eligible ) continue;

        M.set_flag( 2, s, SimplexFlag::SimplexFlagDirichlet );
        
        M.set_flag( 1, M.get_subsimplex( 2, 1, s, 0 ), SimplexFlag::SimplexFlagDirichlet );
        M.set_flag( 1, M.get_subsimplex( 2, 1, s, 1 ), SimplexFlag::SimplexFlagDirichlet );
        M.set_flag( 1, M.get_subsimplex( 2, 1, s, 2 ), SimplexFlag::SimplexFlagDirichlet );
        
        M.set_flag( 0, M.get_subsimplex( 2, 0, s, 0 ), SimplexFlag::SimplexFlagDirichlet );
        M.set_flag( 0, M.get_subsimplex( 2, 0, s, 1 ), SimplexFlag::SimplexFlagDirichlet );
        M.set_flag( 0, M.get_subsimplex( 2, 0, s, 2 ), SimplexFlag::SimplexFlagDirichlet );
        
    }
    
    M.check_dirichlet_flags(false);
    
    LOG << "Prepare scalar fields for testing..." << nl;


    

    
    
    // std::cout << M.outputTikZ(true);
    

    

    LOG << "Unit Test: 3D transformed Poisson Problem (Lagrange matrices)" << nl;

    const int min_l = 1; 
    const int max_l = 4;

    std::vector<MeshSimplicial3D>    meshes;
    std::vector<FloatVector>         solutions;
    std::vector<ConvergenceTable>    contables;
    
    ConvergenceTable contable_pastone("Mass errror, one before");
    ConvergenceTable contable_pasttwo("Mass errror, two before");   

    contable_pastone << "level" << "u_error" << "du_error" << nl;
    contable_pasttwo << "level" << "u_error" << "du_error" << nl;


    for( int l = 0; l < min_l; l++ )
        M.uniformrefinement();

    for( int l = min_l; l <= max_l; l++ ){
        
        display_mallinfo();

        LOG << "Level: " << l << "/" << max_l << nl;
        LOG << "# T/F/E/V: " << M.count_tetrahedra() << "/" << M.count_faces() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        {
            
            int r = 1;
            int w = r;

            LOGPRINTF("Polynomial degrees: r=%d w=%d \n", r, w );
            
            LOG << "...assemble scalar mass matrices" << nl;
    
            auto mass = LagrangeCoefficientMassMatrix( M, r, w, weight_scalar );
            // auto mass = MatrixCSR( LagrangeMassMatrix( M, r ) );

            LOG << "...assemble vector mass matrix" << nl;
    
            auto stiffness = LagrangeCoefficientStiffnessMatrix( M, r, w, weight_vector );
            // auto stiffness = MatrixCSR( LagrangeStiffnessMatrix( M, r ) );
            
            // display_mallinfo();

            mass.compressentries();
            stiffness.compressentries();

            LOG << "Size of mass matrix: " << sizeof(int) * ( mass.getdimout() + mass.getnumberofentries() ) + sizeof(Float) * mass.getnumberofentries() << nl;
            LOG << "Size of mass matrix: " << sizeof(int) * ( stiffness.getdimout() + stiffness.getnumberofentries() ) + sizeof(Float) * stiffness.getnumberofentries() << nl;

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
                for( int s = 0; s < M.count_tetrahedra(); s++ )
                {
                    auto measure = M.getMeasure(3,s);

                    auto midpoint = M.get_midpoint(3,s);

                    auto f = parametric_f(midpoint);

                    auto contribution = f * measure / 4.;

                    for( int vi = 0; vi < 4; vi++ )
                    {
                        int v = M.get_subsimplex( 3, 0, s, vi );
                        
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

                LOG << "Size of stiffness matrix: " << sizeof(int) * ( stiffness.getdimout() + stiffness.getnumberofentries() ) + sizeof(Float) * stiffness.getnumberofentries() << nl;
                
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

                
                if( true and r == 1 )
                {
                    fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
                    VTKWriter vtk( M, fs, getbasename(__FILE__) );
                    
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
