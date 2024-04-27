

/**/

#include <ostream>
// #include <fstream>
// #include <iomanip>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../sparse/rainbow.hpp"
// #include "../../solver/chebyshev.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Solve SPD system: CSR Solvers" );

int main( int argc, char *argv[] )
{
        LOG << "Unit Test: " << TestName << nl;
        
        // LOG << std::setprecision(5);

        if(true){

            ConvergenceTable contable_sol("L2 Error");
            ConvergenceTable contable_res("L2 Residual");
            ConvergenceTable contable_num("Iteration percentage");
            ConvergenceTable contable_sec("Time");
            
            contable_sol.print_rowwise_instead_of_columnwise = true;
            contable_res.print_rowwise_instead_of_columnwise = true;
            contable_num.print_rowwise_instead_of_columnwise = true;
            contable_sec.print_rowwise_instead_of_columnwise = true;
            
            bool do_cgmpp      = false;
            bool do_crmpp_expl = false;
            bool do_crmpp_robt = false;
            bool do_crmpp_fast = false;
            bool do_minres     = false;
            bool do_herzog     = false;
            //
            bool do_cgm_csr                = false;
            bool do_crm_csr                = false;
            bool do_crm_csrtextbook        = false;
            bool do_minres_csr             = false;
            bool do_whatever_csr           = false;
            bool do_cgm_diagonal_csr       = false;
            bool do_cgm_ssor_csr           = false;
            bool do_chebyshev_diagonal_csr = false;
            bool do_cgm_rainbow_csr        = false;

            do_cgmpp      = true;
            do_crmpp_expl = true;
            do_crmpp_robt = true;
            do_crmpp_fast = true;
            do_minres     = true;
            do_herzog     = true;
            //
            do_cgm_csr                = true;
            do_crm_csr                = true;
            do_crm_csrtextbook        = true;
            do_minres_csr             = true;
            do_whatever_csr           = true;
            do_cgm_diagonal_csr       = true;
            do_cgm_ssor_csr           = true;
            do_chebyshev_diagonal_csr = true;
            do_cgm_rainbow_csr        = true;

            
            
            
            contable_sol << "Index";
            if( do_cgmpp      ) contable_sol << "CGM++"      ;
            if( do_crmpp_expl ) contable_sol << "CRM++(expl)";
            if( do_crmpp_robt ) contable_sol << "CRM++(robt)";
            if( do_crmpp_fast ) contable_sol << "CRM++(fast)";
            if( do_minres     ) contable_sol << "MINRES"     ;
            if( do_herzog     ) contable_sol << "HERZOG"     ;
            //
            if( do_cgm_csr )                contable_sol << "CGMcsr"       ;
            if( do_crm_csr )                contable_sol << "CRMcsr"       ;
            if( do_crm_csrtextbook )        contable_sol << "CRMcsr_tb"    ;
            if( do_minres_csr )             contable_sol << "MINREScsr"    ;
            if( do_whatever_csr )           contable_sol << "WHATEVER"     ;
            if( do_cgm_diagonal_csr )       contable_sol << "CGMcsr_diag"  ;
            if( do_cgm_ssor_csr )           contable_sol << "CGMcsr_ssor"  ;
            if( do_chebyshev_diagonal_csr ) contable_sol << "Chebyshev_csr";
            if( do_cgm_rainbow_csr )        contable_sol << "CGMcsr_rainbow";
            contable_sol << nl;
            
            contable_res << "Index";
            if( do_cgmpp      ) contable_res << "CGM++"      ;
            if( do_crmpp_expl ) contable_res << "CRM++(expl)";
            if( do_crmpp_robt ) contable_res << "CRM++(robt)";
            if( do_crmpp_fast ) contable_res << "CRM++(fast)";
            if( do_minres     ) contable_res << "MINRES"     ;
            if( do_herzog     ) contable_res << "HERZOG"     ;
            //
            if( do_cgm_csr )                contable_res << "CGMcsr"       ;
            if( do_crm_csr )                contable_res << "CRMcsr"       ;
            if( do_crm_csrtextbook )        contable_res << "CRMcsr_tb"    ;
            if( do_minres_csr )             contable_res << "MINREScsr"    ;
            if( do_whatever_csr )           contable_res << "WHATEVER"     ;
            if( do_cgm_diagonal_csr )       contable_res << "CGMcsr_diag"  ;
            if( do_cgm_ssor_csr )           contable_res << "CGMcsr_ssor"  ;
            if( do_chebyshev_diagonal_csr ) contable_res << "Chebyshev_csr";
            if( do_cgm_rainbow_csr )        contable_res << "CGMcsr_rainbow";
            contable_res << nl;

            contable_num << "Index";
            if( do_cgmpp      ) contable_num << "CGM++"      ;
            if( do_crmpp_expl ) contable_num << "CRM++(expl)";
            if( do_crmpp_robt ) contable_num << "CRM++(robt)";
            if( do_crmpp_fast ) contable_num << "CRM++(fast)";
            if( do_minres     ) contable_num << "MINRES"     ;
            if( do_herzog     ) contable_num << "HERZOG"     ;
            //
            if( do_cgm_csr )                contable_num << "CGMcsr"       ;
            if( do_crm_csr )                contable_num << "CRMcsr"       ;
            if( do_crm_csrtextbook )        contable_num << "CRMcsr_tb"    ;
            if( do_minres_csr )             contable_num << "MINREScsr"    ;
            if( do_whatever_csr )           contable_num << "WHATEVER"     ;
            if( do_cgm_diagonal_csr )       contable_num << "CGMcsr_diag"  ;
            if( do_cgm_ssor_csr )           contable_num << "CGMcsr_ssor"  ;
            if( do_chebyshev_diagonal_csr ) contable_num << "Chebyshev_csr";
            if( do_cgm_rainbow_csr )        contable_num << "CGMcsr_rainbow";
            contable_num << nl;


            contable_sec << "Index";
            if( do_cgmpp      ) contable_sec << "CGM++"      ;
            if( do_crmpp_expl ) contable_sec << "CRM++(expl)";
            if( do_crmpp_robt ) contable_sec << "CRM++(robt)";
            if( do_crmpp_fast ) contable_sec << "CRM++(fast)";
            if( do_minres     ) contable_sec << "MINRES"     ;
            if( do_herzog     ) contable_sec << "HERZOG"     ;
            //
            if( do_cgm_csr )                contable_sec << "CGMcsr"       ;
            if( do_crm_csr )                contable_sec << "CRMcsr"       ;
            if( do_crm_csrtextbook )        contable_sec << "CRMcsr_tb"    ;
            if( do_minres_csr )             contable_sec << "MINREScsr"    ;
            if( do_whatever_csr )           contable_sec << "WHATEVER"     ;
            if( do_cgm_diagonal_csr )       contable_sec << "CGMcsr_diag"  ;
            if( do_cgm_ssor_csr )           contable_sec << "CGMcsr_ssor"  ;
            if( do_chebyshev_diagonal_csr ) contable_sec << "Chebyshev_csr";
            if( do_cgm_rainbow_csr )        contable_sec << "CGMcsr_rainbow";
            contable_sec << nl;

            

            const std::vector<int> Ns = { 4, 8, 16, 32 };

            for( const int N : Ns ){
                
                LOG << "Level: " << N << nl;
                
                {
                    
                    LOG << "...assemble matrix" << nl;

                    std::vector< SparseMatrix::MatrixEntry > entries;

                    const Float h = Float(1.) / (N+1);
                    const Float h2 = h * h;

                    for( int e = 0; e < N*N; e++ )
                    {
                        int x = e / N;
                        int y = e % N;
                        assert( e == x * N + y );

                        entries.push_back({ x * N + y, x * N + y, Float(4.) / h2 });

                        if( x != 0   ) entries.push_back({ x * N + y, (x-1) * N + y,   Float(-1.) / h2 });
                        if( x != N-1 ) entries.push_back({ x * N + y, (x+1) * N + y,   Float(-1.) / h2 });
                        if( y != 0   ) entries.push_back({ x * N + y, (x  ) * N + y-1, Float(-1.) / h2 });
                        if( y != N-1 ) entries.push_back({ x * N + y, (x  ) * N + y+1, Float(-1.) / h2 });
                        
                    }

                    auto system_prelim = SparseMatrix( N*N, N*N, entries );
                    system_prelim.sortentries();
                    auto system = MatrixCSR( system_prelim );


                    LOG << "...create solutions and right-hand sides" << nl;

                    const int T = 5;

                    std::vector<FloatVector> sols;
                    std::vector<FloatVector> rhss;

                    for( int t = 0; t < T; t++ )
                    {
                        FloatVector sol( N * N );
                        sol.random();
                        sol.normalize();
                        FloatVector rhs = system * sol;
                        sols.push_back( sol );
                        rhss.push_back( rhs );
                    }



                    for( int t = 0; t < T; t++ )
                    {

                        const auto& sol = sols[t];
                        const auto& rhs = rhss[t];
                        
                        contable_sol << static_cast<Float>(N);
                        contable_res << static_cast<Float>(N);
                        contable_num << static_cast<Float>(N);
                        contable_sec << static_cast<Float>(N);

                        if( do_cgmpp )
                        {
                            LOG << "CGM C++" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateGradientMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve( mysol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }

                        if( do_crmpp_expl )
                        {
                            LOG << "CRM C++ (explicit)" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = 4 * mysol.getdimension();
                            
                            timestamp start = timestampnow();
                            Solver.solve_explicit( mysol, rhs );
                            timestamp end = timestampnow();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }

                        if( do_crmpp_robt )
                        {
                            LOG << "CRM C++ (robust)" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = 4 * mysol.getdimension();
                            
                            timestamp start = timestampnow();
                            Solver.solve_robust( mysol, rhs );
                            timestamp end = timestampnow();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }

                        if( do_crmpp_fast )
                        {
                            LOG << "CRM C++ (fast)" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = 4 * mysol.getdimension();
                            
                            timestamp start = timestampnow();
                            Solver.solve_fast( mysol, rhs );
                            timestamp end = timestampnow();
                            
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }

                        if( do_minres )
                        {
                            LOG << "MINRES C++" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            MinimumResidualMethod Solver( system );
                            Solver.print_modulo        = 4 * mysol.getdimension();
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve( mysol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }

                        if( do_herzog )
                        {
                            LOG << "HERZOG SOODHALTER C++" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            HerzogSoodhalterMethod Solver( system );
                            Solver.print_modulo        = 4 * mysol.getdimension();
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     4 * mysol.getdimension();
                            timestamp start = timestampnow();
                            Solver.solve( mysol, rhs );
                            timestamp end = timestampnow();
                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;

                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( Solver.recent_iteration_count ) / Solver.max_iteration_count;
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }
                        
                        
                        




                        if( do_cgm_csr )
                        {
                            LOG << "CGM - CSR Classic" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            ConjugateGradientSolverCSR( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }

                        if( do_crm_csr )
                        {
                            LOG << "CRM - CSR Classic" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            ConjugateResidualSolverCSR( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }

                        if( do_crm_csrtextbook )
                        {
                            LOG << "CRM - CSR Textbook" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            ConjugateResidualSolverCSR_textbook( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }

                        if( do_minres_csr )
                        {
                            LOG << "MINRES CSR" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            MINRESCSR( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }


                        if( do_whatever_csr )
                        {
                            LOG << "WHATEVER CSR" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            WHATEVER( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }


                        if( do_cgm_diagonal_csr )
                        {
                            LOG << "CGM diagonal preconditioner CSR" << nl;
                        
                            DiagonalOperator invprecon = InverseDiagonalPreconditioner( system_prelim );
//                             invprecon.setentries( 1. );
                            assert( invprecon.getdiagonal().isfinite() );
                            assert( invprecon.getdiagonal().isnonnegative() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            ConjugateGradientSolverCSR_DiagonalPreconditioner( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                invprecon.getdiagonal().raw()
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }
                        
                        
                        if( do_cgm_ssor_csr )
                        {
                            LOG << "CGM SSOR preconditioner CSR" << nl;
                        
                            FloatVector diagonal = system.getDiagonal();
                            assert( diagonal.isfinite() );
                            assert( diagonal.isnonnegative() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            ConjugateGradientSolverCSR_SSOR( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                diagonal.raw(),
                                0.9123456789
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }
                        
                        
                        if( do_chebyshev_diagonal_csr )
                        {
                            LOG << "CHEBYSHEV CSR" << nl;
                        
                            DiagonalOperator invprecon = InverseDiagonalPreconditioner( system_prelim );
                            assert( invprecon.getdiagonal().isfinite() );
                            assert( invprecon.getdiagonal().ispositive() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            ChebyshevIteration_DiagonalPreconditioner( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                invprecon.getdiagonal().raw(),
                                0.,
                                system.eigenvalueupperbound()
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }
                        

                        if( do_cgm_rainbow_csr )
                        {
                            LOG << "CGM Rainbow-SSOR preconditioner CSR" << nl;
                        
                            FloatVector diagonal = system.getDiagonal();
                            assert( diagonal.isfinite() );
                            assert( diagonal.isnonnegative() );
                            
                            FloatVector mysol( N*N );
                            mysol.zero();
                            FloatVector residual( rhs );

                            Rainbow rainbow( system );

                            timestamp start = timestampnow();
                            int recent_iteration_count =
                            ConjugateGradientSolverCSR_Rainbow( 
                                mysol.getdimension(), 
                                mysol.raw(), 
                                rhs.raw(), 
                                system.getA(), system.getC(), system.getV(),
                                residual.raw(),
                                desired_precision,
                                0,
                                diagonal.raw(),
                                0.9123456789,
                                rainbow.num_colors, rainbow.F.data(), rainbow.B.data(), rainbow.R.data()
                            );
                            timestamp end = timestampnow();

                            LOG << "\t\t\t Time: " << timestamp2measurement( end - start ) << nl;
                            
                            auto stat_sol = Float( ( sol - mysol ).norm() );
                            auto stat_res = Float( ( system * mysol - rhs ).norm() );
                            auto stat_num = Float( recent_iteration_count )/ (N*N); 
                            contable_sol << stat_sol;
                            contable_res << stat_res;
                            contable_num << stat_num;
                            contable_sec << Float( end - start );
                        }
                        
                        
                        contable_sol << nl;
                        contable_res << nl;
                        contable_num << nl;
                        contable_sec << nl;
                        
                    }
                    
                    contable_sol.lg( false );
                    contable_res.lg( false );
                    contable_num.lg( false );
                    contable_sec.lg( false );

                    }

                
                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test: " << TestName << nl;
        
        return 0;
}
