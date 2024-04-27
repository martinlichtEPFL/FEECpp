

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
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../solver/inv.hpp"
#include "../../solver/systemsparsesolver.hpp"
#include "../../solver/nullspace.hpp"
// #include "../../solver/cgm.hpp"
// #include "../../solver/crm.hpp"
// #include "../../solver/pcrm.hpp"
// #include "../../solver/minres.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/global.interpol.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

const Float mass_threshold_for_small_vectors = 1e-6;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: Nullspace computation (2D) Hodge-Laplacian" << nl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D Mx = StandardSquare2D_tiles3x3();
            
            Mx.automatic_dirichlet_flags();

            // LOG << Mx << nl;

            // Seitenmitten: 2, 11, 19, 31
            Mx.set_flag( 1,  2, SimplexFlag::SimplexFlagNull );
            Mx.set_flag( 1, 11, SimplexFlag::SimplexFlagNull );
            Mx.set_flag( 1, 19, SimplexFlag::SimplexFlagNull );
            Mx.set_flag( 1, 31, SimplexFlag::SimplexFlagNull );

            // Links: 1, 11, 21 Rechts: 9, 19, 29
            // Mx.set_flag( 1,  1, SimplexFlag::SimplexFlagNull );
            // Mx.set_flag( 1, 11, SimplexFlag::SimplexFlagNull );
            // Mx.set_flag( 1, 21, SimplexFlag::SimplexFlagNull );
            // Mx.set_flag( 1,  9, SimplexFlag::SimplexFlagNull );
            // Mx.set_flag( 1, 19, SimplexFlag::SimplexFlagNull );
            // Mx.set_flag( 1, 29, SimplexFlag::SimplexFlagNull );
            
            // Mx.set_flag( 0, 4, SimplexFlag::SimplexFlagNull );
            // Mx.set_flag( 0, 8, SimplexFlag::SimplexFlagNull );
            // Mx.set_flag( 0, 7, SimplexFlag::SimplexFlagNull );
            // Mx.set_flag( 0,11, SimplexFlag::SimplexFlagNull );
            



            Mx.check();
            
            
            
            MeshSimplicial2D M;
            
            for( int i = 0; i < 1; i++ )
            {
                auto M2 = Mx;
                M2.getcoordinates().shift( FloatVector{ i * 3.0, 0.0 } );
                M.merge( M2 );
            }
                        
            
            std::function<FloatVector(const FloatVector&)> constant_one
                = [](const FloatVector& vec) -> FloatVector{
                        assert( vec.getdimension() == 2 );
                        return FloatVector({ 1. });
                    };
            

            LOG << "Nullspace computation" << nl;

            ConvergenceTable contable("Number of nullvectors");
            
            contable << "#nullvec" << nl;

            const Float desired_precision = 100 * machine_epsilon;
            

            const int min_l = 0; 
            
            const int max_l = 5;
            
            const int min_r = 1; 
            
            const int max_r = 1;
            
            const int max_number_of_candidates = 4;

            const int max_number_of_purifications = 1;

            assert( 0 <= min_l and min_l <= max_l );
            assert( 0 <= min_r and min_r <= max_r );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ )
            {
                
                LOG << "Level: " << l << "/" << max_l << nl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                for( int r = min_r; r <= max_r; r++ )
                {
                    
                    LOG << "Polynomial degree: " << r << "/" << max_r << nl;
                    
                    LOG << "...assemble partial matrices" << nl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r   );
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix volume_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 );

                    SparseMatrix vector_elevationmatrix   = FEECBrokenElevationMatrix( M, M.getinnerdimension(), 1, r-1, 1);
                    SparseMatrix vector_elevationmatrix_t = vector_elevationmatrix.getTranspose();

                    SparseMatrix scalar_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r   );
                    SparseMatrix scalar_incmatrix_t = scalar_incmatrix.getTranspose();

                    SparseMatrix vector_incmatrix   = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r   );
                    SparseMatrix vector_incmatrix_t = vector_incmatrix.getTranspose();

                    SparseMatrix scalar_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );
                    SparseMatrix scalar_diffmatrix_t = scalar_diffmatrix.getTranspose();

                    SparseMatrix vector_diffmatrix   = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r );
                    SparseMatrix vector_diffmatrix_t = vector_diffmatrix.getTranspose();


                    LOG << "... full matrices" << nl;
            
                    auto mass = vector_incmatrix_t * vector_massmatrix * vector_incmatrix;

                    auto mat_A  = scalar_incmatrix_t & scalar_massmatrix & scalar_incmatrix;
                    mat_A.sortandcompressentries();
                    
                    LOG << vector_elevationmatrix.getdimin() << space << vector_elevationmatrix.getdimout() << nl;

                    auto mat_Bt = scalar_incmatrix_t & scalar_diffmatrix_t & vector_elevationmatrix_t & vector_massmatrix & vector_incmatrix; // upper right
                    mat_Bt.sortandcompressentries();
                    
                    auto mat_B = mat_Bt.getTranspose(); //volume_incmatrix_t & volume_massmatrix & diffmatrix & vector_incmatrix; // lower bottom
                    mat_B.sortandcompressentries();
                    
                    auto mat_C  = vector_incmatrix_t & vector_diffmatrix_t & volume_massmatrix & vector_diffmatrix & vector_incmatrix;
                    mat_C.sortandcompressentries();
                    
                    auto A  = MatrixCSR( mat_A  );
                    auto Bt = MatrixCSR( mat_Bt );
                    auto B  = MatrixCSR( mat_B  );
                    auto C  = MatrixCSR( mat_C  );
                    
                    auto Z  = MatrixCSR( mat_B.getdimout(), mat_B.getdimout() ); // zero matrix
                    
                    auto SystemMatrix = C + B * inv(A,100*machine_epsilon,-1) * Bt;
                    
                    
                    
                    
                    
                    
                    
                    std::vector<FloatVector> nullvectorgallery;
                    
                    
                    for( int no_candidate = 0; no_candidate < max_number_of_candidates; no_candidate++ )
                    {
                        
                        FloatVector candidate( Bt.getdimin(), 0. ); 
                        candidate.random(); 
                        candidate.normalize(mass);
                        
                        {
                            for( int s = 0; s < 2; s++ )
                            for( const auto& nullvector : nullvectorgallery ) {
                                Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                                candidate = candidate - alpha * nullvector;
                            }
                            
                            Float reduced_mass = candidate.norm(mass);
                            LOG << "\t\t\t Preprocessed mass: " << reduced_mass << nl;
                            
                            if( reduced_mass < mass_threshold_for_small_vectors ) {
                                LOG << "**** The candidate already has very small mass" << nl;
//                                 continue;
                            }
                        }
                        
                        
                        /* reduce the candidate to its nullspace component */
                        {
                            const FloatVector rhs( Bt.getdimin(), 0. );
                        
                            FloatVector residual( rhs );
                            
                            for( int t = 0; t < max_number_of_purifications; t++ )
                            {
                                
                                auto& X = SystemMatrix;

                                HodgeConjugateResidualSolverCSR_SSOR(
                                    B.getdimout(), 
                                    A.getdimout(), 
                                    candidate.raw(), 
                                    rhs.raw(), 
                                    A.getA(),   A.getC(),  A.getV(), 
                                    B.getA(),   B.getC(),  B.getV(), 
                                    Bt.getA(), Bt.getC(), Bt.getV(), 
                                    C.getA(),   C.getC(),  C.getV(), 
                                    residual.raw(),
                                    desired_precision,
                                    0,
                                    100*machine_epsilon,
                                    -1
                                );
                                
                                // ConjugateResidualSolverCSR( 
                                //     candidate.getdimension(), 
                                //     candidate.raw(), 
                                //     rhs.raw(), 
                                //     C.getA(), C.getC(), C.getV(),
                                //     residual.raw(),
                                //     desired_precision,
                                //     0
                                // );
                                
                                // ConjugateResidualMethod Solver( X );
                                // Solver.tolerance           = desired_precision;
                                // Solver.print_modulo        = 0;
                                // Solver.max_iteration_count = 1 * candidate.getdimension();
                                // Solver.solve( candidate, rhs );
                            
                                assert( candidate.isfinite() );
                                
                                LOG << "\t\t\t (eucl) delta:     " << ( residual - rhs + X * candidate ).norm() << nl;
                                LOG << "\t\t\t (mass) delta:     " << ( residual - rhs + X * candidate ).norm( mass ) << nl;
                                LOG << "\t\t\t (eucl) res:       " << residual.norm() << nl;
                                LOG << "\t\t\t (mass) res:       " << residual.norm( mass ) << nl;
                                LOG << "\t\t\t (eucl) x:         " << candidate.norm() << nl;
                                LOG << "\t\t\t (mass) x:         " << candidate.norm( mass ) << nl;
                                LOG << "\t\t\t (eucl) Ax:        " << ( X * candidate ).norm() << nl;
                                LOG << "\t\t\t (mass) Ax:        " << ( X * candidate ).norm( mass ) << nl;
                                LOG << "\t\t\t (eucl) b - Ax:    " << ( X * candidate - rhs ).norm() << nl;
                                LOG << "\t\t\t (mass) b - Ax:    " << ( X * candidate - rhs ).norm( mass ) << nl;
                                
                                candidate.normalize( mass );
                                
                                assert( candidate.isfinite() );
                                
                                LOG << "\t\t\t (norm eucl) x:         " << candidate.norm() << nl;
                                LOG << "\t\t\t (norm mass) x:         " << candidate.norm( mass ) << nl;
                                LOG << "\t\t\t (norm eucl) Ax:        " << ( X * candidate ).norm() << nl;
                                LOG << "\t\t\t (norm mass) Ax:        " << ( X * candidate ).norm( mass ) << nl;

                            }
                        }
                        
                        
                        /* Is it nullspace at all? */
                        
                        candidate.normalize(mass);
                        
                        Float residual_mass = ( SystemMatrix * candidate ).norm(mass);
                        
                        LOG << "\t\t\t Numerical residual (after normalizing): " << residual_mass << nl;
                        
                        if( residual_mass > mass_threshold_for_small_vectors ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because not nullspace enough!" << nl;
                            continue;
                        }
                        
                        /* Gram-Schmidt */
                        
                        for( int s = 0; s < 2; s++ )
                        for( const auto& nullvector : nullvectorgallery ) {
                            Float alpha = (mass*candidate*nullvector) / (mass*nullvector*nullvector);
                            candidate = candidate - alpha * nullvector;
                        }
                        
                        Float reduced_mass = candidate.norm(mass);
                        LOG << "\t\t\t Reduced mass: " << reduced_mass << nl;
                        LOG << "\t\t\t Numerical residual (after Gram-Schmidt): " << ( SystemMatrix * candidate ).norm(mass) << nl;
                        
                        if( reduced_mass < mass_threshold_for_small_vectors ) {
                            LOG << "!!!!!!!!!!!!!Discard vector because mass is too small!" << nl;
                            continue;
                        }
                        
                        candidate.normalize(mass);

                        assert( candidate.isfinite() );
                        
                        LOG << "Accept vector #" << nullvectorgallery.size() + 1 << nl;
                    
                        
                        nullvectorgallery.push_back( candidate );
                    }
                    
                    
                    
                    LOG << "How much nullspace are our vectors?" << nl;
                    for( const auto& nullvector : nullvectorgallery ) {
                        Float mass_norm = ( SystemMatrix * nullvector ).norm(mass);
                        Assert( mass_norm < mass_threshold_for_small_vectors, mass_norm, mass_threshold_for_small_vectors );
                        // LOGPRINTF( "% 10.5Le\t", (long double)mass_norm );
                        LOG << mass_norm << tab;
                    }
                    LOG << nl;
                    
                    LOG << "How orthonormal are our vectors?" << nl;
                    for( int n1 = 0; n1 < nullvectorgallery.size(); n1++ ) {
                        for( int n2 = 0; n2 < nullvectorgallery.size(); n2++ ) {
                            auto nullvector1 = nullvectorgallery[n1];
                            auto nullvector2 = nullvectorgallery[n2];
                            Float mass_prod = mass * nullvector1 * nullvector2;
                            // LOGPRINTF( "% 10.5Le\t", (long double)mass_prod );
                            LOG << mass_prod << tab;
                            if( n1 != n2 ) assert( is_numerically_small( mass_prod ) );
                            
                        }
                        LOG << nl;
                    }
                    
                    
                    
                    contable << static_cast<Float>(nullvectorgallery.size());  

                    
                    const auto interpol_matrix = FEECBrokenInterpolationMatrix( M, M.getinnerdimension(), 1, 0, r );

                    for( const auto& nullvector : nullvectorgallery )
                    {
                
                        fstream fs( experimentfile(getbasename(__FILE__)), std::fstream::out );
            
                        VTKWriter vtk( M, fs, getbasename(__FILE__) );

                        auto reduced_nullvector = interpol_matrix * vector_incmatrix * nullvector;

                        vtk.writeCellVectorData_barycentricgradients( reduced_nullvector, "nullvector_Hcurl" , 1.0 );
                        
                        fs.close();
                
                    } 

                }

                if( l != max_l ) { LOG << "Refinement..." << nl; M.uniformrefinement(); }

                contable << nl;
                
                contable.lg();
        
            } 
            
            contable.lg();
        
        }
        
        
        
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}



