
#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/utilities.hpp"


int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: Compare numerical solvers CRM vs MINRES\n           for Solution of Dirichlet Problem" << nl;
        
        // LOG << std::setprecision(10);

        if(true){

            LOG << "Initial mesh..." << nl;
            
            MeshSimplicial2D M = StandardSquare2D();
            
            M.check();
            
            M.automatic_dirichlet_flags();
            
            M.check_dirichlet_flags();
            
            LOG << "Prepare scalar fields for testing..." << nl;
            
            
            const Float xfeq = 1.;
            const Float yfeq = 1.;

            std::function<FloatVector(const FloatVector&)> experiment_rhs = 
                [=](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 2 );
                    return FloatVector({ 
                        
                        xfeq*xfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                        +
                        yfeq*yfeq * Constants::fourpisquare * std::sin( xfeq * Constants::twopi * vec[0] ) * std::sin( yfeq * Constants::twopi * vec[1] )
                     });
                };
            

            
            
            

            // ConvergenceTable contable_sol("L2 Error");
            ConvergenceTable contable_res("L2 Residual");
            ConvergenceTable contable_num("Iteration percentage");
            ConvergenceTable contable_sec("Time");

            // contable_sol.print_rowwise_instead_of_columnwise = true;
            contable_res.print_rowwise_instead_of_columnwise = true;
            contable_num.print_rowwise_instead_of_columnwise = true;
            contable_sec.print_rowwise_instead_of_columnwise = true;
            
            bool do_cgmpp      = true;
            bool do_crmpp_expl = true;
            bool do_crmpp_robt = true;
            bool do_crmpp_fast = true;
            bool do_minres     = true;
            bool do_herzog     = true;
            //
            // bool do_cgm_csr                = true;
            // bool do_crm_csr                = true;
            // bool do_crm_csrtextbook        = true;
            // bool do_minres_csr             = true;
            // bool do_whatever_csr           = false;
            // bool do_cgm_diagonal_csr       = true;
            // bool do_cgm_ssor_csr           = true;
            // bool do_chebyshev_diagonal_csr = true;

            // if( do_cgmpp      ) contable_sol << "CGM++"      ;
            // if( do_crmpp_expl ) contable_sol << "CRM++(expl)";
            // if( do_crmpp_robt ) contable_sol << "CRM++(robt)";
            // if( do_crmpp_fast ) contable_sol << "CRM++(fast)";
            // if( do_minres     ) contable_sol << "MINRES"     ;
            // if( do_herzog     ) contable_sol << "HERZOG"     ;
            // //
            // if( do_cgm_csr )                contable_sol << "CGMcsr"       ;
            // if( do_crm_csr )                contable_sol << "CRMcsr"       ;
            // if( do_crm_csrtextbook )        contable_sol << "CRMcsr_tb"    ;
            // if( do_minres_csr )             contable_sol << "MINREScsr"    ;
            // if( do_whatever_csr )           contable_sol << "WHATEVER"     ;
            // if( do_cgm_diagonal_csr )       contable_sol << "CGMcsr_diag"  ;
            // if( do_cgm_ssor_csr )           contable_sol << "CGMcsr_ssor"  ;
            // if( do_chebyshev_diagonal_csr ) contable_sol << "Chebyshev_csr";
            // contable_sol << nl;
            
            if( do_cgmpp      ) contable_res << "CGM++"      ;
            if( do_crmpp_expl ) contable_res << "CRM++(expl)";
            if( do_crmpp_robt ) contable_res << "CRM++(robt)";
            if( do_crmpp_fast ) contable_res << "CRM++(fast)";
            if( do_minres     ) contable_res << "MINRES"     ;
            if( do_herzog     ) contable_res << "HERZOG"     ;
            //
            // if( do_cgm_csr )                contable_res << "CGMcsr"       ;
            // if( do_crm_csr )                contable_res << "CRMcsr"       ;
            // if( do_crm_csrtextbook )        contable_res << "CRMcsr_tb"    ;
            // if( do_minres_csr )             contable_res << "MINREScsr"    ;
            // if( do_whatever_csr )           contable_res << "WHATEVER"     ;
            // if( do_cgm_diagonal_csr )       contable_res << "CGMcsr_diag"  ;
            // if( do_cgm_ssor_csr )           contable_res << "CGMcsr_ssor"  ;
            // if( do_chebyshev_diagonal_csr ) contable_res << "Chebyshev_csr";
            contable_res << nl;

            if( do_cgmpp      ) contable_num << "CGM++"      ;
            if( do_crmpp_expl ) contable_num << "CRM++(expl)";
            if( do_crmpp_robt ) contable_num << "CRM++(robt)";
            if( do_crmpp_fast ) contable_num << "CRM++(fast)";
            if( do_minres     ) contable_num << "MINRES"     ;
            if( do_herzog     ) contable_num << "HERZOG"     ;
            //
            // if( do_cgm_csr )                contable_num << "CGMcsr"       ;
            // if( do_crm_csr )                contable_num << "CRMcsr"       ;
            // if( do_crm_csrtextbook )        contable_num << "CRMcsr_tb"    ;
            // if( do_minres_csr )             contable_num << "MINREScsr"    ;
            // if( do_whatever_csr )           contable_num << "WHATEVER"     ;
            // if( do_cgm_diagonal_csr )       contable_num << "CGMcsr_diag"  ;
            // if( do_cgm_ssor_csr )           contable_num << "CGMcsr_ssor"  ;
            // if( do_chebyshev_diagonal_csr ) contable_num << "Chebyshev_csr";
            contable_num << nl;
            
            if( do_cgmpp      ) contable_sec << "CGM++"      ;
            if( do_crmpp_expl ) contable_sec << "CRM++(expl)";
            if( do_crmpp_robt ) contable_sec << "CRM++(robt)";
            if( do_crmpp_fast ) contable_sec << "CRM++(fast)";
            if( do_minres     ) contable_sec << "MINRES"     ;
            if( do_herzog     ) contable_sec << "HERZOG"     ;
            //
            // if( do_cgm_csr )                contable_sec << "CGMcsr"       ;
            // if( do_crm_csr )                contable_sec << "CRMcsr"       ;
            // if( do_crm_csrtextbook )        contable_sec << "CRMcsr_tb"    ;
            // if( do_minres_csr )             contable_sec << "MINREScsr"    ;
            // if( do_whatever_csr )           contable_sec << "WHATEVER"     ;
            // if( do_cgm_diagonal_csr )       contable_sec << "CGMcsr_diag"  ;
            // if( do_cgm_ssor_csr )           contable_sec << "CGMcsr_ssor"  ;
            // if( do_chebyshev_diagonal_csr ) contable_sec << "Chebyshev_csr";
            contable_sec << nl;
            

            

            const int min_l = 0; 
            const int max_l = 6;

            assert( 0 <= min_l and min_l <= max_l );
            
            for( int l = 0; l < min_l; l++ )
                M.uniformrefinement();

            for( int l = min_l; l <= max_l; l++ ){
                
                LOG << "Level: " << l << "/" << max_l << nl;
                LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
                
                const int r = 1;
                
                {
                    
                    LOG << "...assemble scalar mass matrix" << nl;
            
                    SparseMatrix scalar_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                    
                    LOG << "...assemble vector mass matrix" << nl;
            
                    SparseMatrix vector_massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 );

                    LOG << "...sort and compress mass matrices" << nl;
            
                    scalar_massmatrix.sortandcompressentries();
                    vector_massmatrix.sortandcompressentries();

                    LOG << "...assemble differential matrix and transpose" << nl;

                    SparseMatrix diffmatrix = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix diffmatrix_t = diffmatrix.getTranspose();

                    LOG << "...assemble inclusion matrix and transpose" << nl;
            
                    SparseMatrix incmatrix = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );

                    SparseMatrix incmatrix_t = incmatrix.getTranspose();

                    LOG << "...compose stiffness and mass matrices" << nl;
            
                    const auto composed_stiffness = incmatrix_t * diffmatrix_t * vector_massmatrix * diffmatrix * incmatrix;
                    const auto composed_mass      = incmatrix_t * scalar_massmatrix * incmatrix;

                    auto opr  = diffmatrix & incmatrix;
                    auto opl  = opr.getTranspose(); 
                    auto stiffness_csr_prelim = opl & ( vector_massmatrix & opr );
                    
                    LOG << "...convert to CSR" << nl;
            
                    stiffness_csr_prelim.sortentries();
                    auto stiffness_csr = MatrixCSR( stiffness_csr_prelim );

                    const auto& stiffness = composed_stiffness;
                    const auto& mass      = composed_mass;
                    
                    LOG << "...matrices done" << nl;

                    {

                        FloatVector sol_original( M.count_simplices(0), 0. );
                        sol_original.random();
                        sol_original.normalize( mass );
                        
                        const auto& function_rhs  = experiment_rhs;
                        FloatVector interpol_rhs  = Interpolation( M, M.getinnerdimension(), 0, r,   function_rhs  );
                        FloatVector rhs_original = incmatrix_t * ( scalar_massmatrix * interpol_rhs );
                        
                        // rhs_original.zero();
                        
                        // const Float desired_precision = sqrt( machine_epsilon );

                        if( do_cgmpp )
                        {
                            LOG << "CGM C++" << nl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            ConjugateGradientMethod Solver( stiffness );
                            Solver.print_modulo        = 0;
                            Solver.tolerance        = desired_precision;
                            Solver.max_iteration_count =     1 * sol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve( sol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_crmpp_expl )
                        {
                            LOG << "CRM C++" << nl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            ConjugateResidualMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.tolerance        = desired_precision;
                            Solver.max_iteration_count =     1 * sol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve_explicit( sol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_crmpp_robt )
                        {
                            LOG << "CRM C++" << nl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            ConjugateResidualMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.tolerance        = desired_precision;
                            Solver.max_iteration_count =     1 * sol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve_robust( sol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_crmpp_fast )
                        {
                            LOG << "CRM C++" << nl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            ConjugateResidualMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.tolerance        = desired_precision;
                            Solver.max_iteration_count =     1 * sol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve_fast( sol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_minres )
                        {
                            LOG << "MINRES C++" << nl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            MinimumResidualMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.tolerance        = desired_precision;
                            Solver.max_iteration_count =     1 * sol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve( sol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }

                        if( do_herzog )
                        {
                            LOG << "HERZOG SOODHALTER C++" << nl;
                        
                            FloatVector sol = sol_original;
                            const FloatVector rhs = rhs_original;
                            HerzogSoodhalterMethod Solver( stiffness );
                            // Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            // Solver.print_modulo        = 1;
                            Solver.print_modulo        = 0;
                            Solver.tolerance        = desired_precision;
                            Solver.max_iteration_count =     1 * sol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve( sol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                            LOG << "Mass of approximate solution: " << sol.norm( mass ) << nl;

                            auto runtime  = static_cast<Float>( end - start );
                            // auto stat_sol = Float( ( sol - ... ).norm() );
                            auto stat_res = Float( ( stiffness * sol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            
                            //contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << runtime;
                        }
                        
                        
                        // contable_sol << nl;
                        contable_res << nl;
                        contable_num << nl;
                        contable_sec << nl;

                        // contable_sol.lg( false );
                        contable_res.lg( false );
                        contable_num.lg( false );
                        contable_sec.lg( false );

                    }
                    
                }

                
                if( l != max_l ){ 
                    LOG << "Refinement..." << nl;
                    M.uniformrefinement();
                    M.shake_interior_vertices();
                }

                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
