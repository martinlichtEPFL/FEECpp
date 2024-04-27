

/**/

#include <ostream>
// #include <fstream>
// #include <iomanip>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../solver/sparsesolver.hpp"
#include "../../solver/iterativesolver.hpp"


using namespace std;

extern const char* TestName;
#define TESTNAME( cstr ) const char* TestName = cstr

TESTNAME( "Solve SPD system: CGM, CRM, MINRES, HerzogSoodhalter" );

int main( int argc, char *argv[] )
{
        LOG << "Unit Test: " << TestName << nl;
        
        // LOG << std::setprecision(5);

        if(true)
        {

            ConvergenceTable contable_sol("L2 Error");
            ConvergenceTable contable_res("L2 Residual");
            ConvergenceTable contable_num("Iteration percentage");

            contable_sol.print_rowwise_instead_of_columnwise = true;
            contable_res.print_rowwise_instead_of_columnwise = true;
            contable_num.print_rowwise_instead_of_columnwise = true;
            
            bool do_cgmpp      = false;
            bool do_crmpp_expl = false;
            bool do_crmpp_robt = false;
            bool do_crmpp_fast = false;
            bool do_minres     = false;
            bool do_herzog     = false;
            
            do_cgmpp      = true;
            do_crmpp_expl = true;
            do_crmpp_robt = true;
            do_crmpp_fast = true;
            do_minres     = true;
            do_herzog     = true;
            
            
            
            
            contable_sol << "Index";
            if( do_cgmpp      ) contable_sol << "CGM++"      ;
            if( do_crmpp_expl ) contable_sol << "CRM++(expl)";
            if( do_crmpp_robt ) contable_sol << "CRM++(robt)";
            if( do_crmpp_fast ) contable_sol << "CRM++(fast)";
            if( do_minres     ) contable_sol << "MINRES"     ;
            if( do_herzog     ) contable_sol << "HERZOG"     ;
            contable_sol << nl;

            contable_res << "Index";
            if( do_cgmpp      ) contable_res << "CGM++"      ;
            if( do_crmpp_expl ) contable_res << "CRM++(expl)";
            if( do_crmpp_robt ) contable_res << "CRM++(robt)";
            if( do_crmpp_fast ) contable_res << "CRM++(fast)";
            if( do_minres     ) contable_res << "MINRES"     ;
            if( do_herzog     ) contable_res << "HERZOG"     ;
            contable_res << nl;

            contable_num << "Index";
            if( do_cgmpp      ) contable_num << "CGM++"      ;
            if( do_crmpp_expl ) contable_num << "CRM++(expl)";
            if( do_crmpp_robt ) contable_num << "CRM++(robt)";
            if( do_crmpp_fast ) contable_num << "CRM++(fast)";
            if( do_minres     ) contable_num << "MINRES"     ;
            if( do_herzog     ) contable_num << "HERZOG"     ;
            contable_num << nl;

            
            const std::vector<int> Ns = { 16, 32 };

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

                    auto system = SparseMatrix( N*N, N*N, entries );


                    LOG << "...create solutions and right-hand sides" << nl;

                    const int T = 10;

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

                        if( do_cgmpp )
                        {
                            LOG << "CGM C++" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateGradientMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count =     mysol.getdimension();
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
                        }

                        if( do_crmpp_expl )
                        {
                            LOG << "CRM C++ (explicit)" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = mysol.getdimension();
                            
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
                        }

                        if( do_crmpp_robt )
                        {
                            LOG << "CRM C++ (robust)" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = mysol.getdimension();
                            
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
                        }

                        if( do_crmpp_fast )
                        {
                            LOG << "CRM C++ (fast)" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            ConjugateResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.max_iteration_count = mysol.getdimension();
                            
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
                        }

                        if( do_minres )
                        {
                            LOG << "MINRES C++" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            MinimumResidualMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     mysol.getdimension();
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
                        }

                        if( do_herzog )
                        {
                            LOG << "HERZOG SOODHALTER C++" << nl;
                        
                            FloatVector mysol( N*N );
                            mysol.zero();
                            HerzogSoodhalterMethod Solver( system );
                            Solver.print_modulo        = mysol.getdimension() / 20;
                            Solver.verbosity        = MinimumResidualMethod::VerbosityLevel::verbose;
                            Solver.max_iteration_count =     mysol.getdimension();
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
                        }
                        
                        
                        contable_sol << nl;
                        contable_res << nl;
                        contable_num << nl;
                        
                    }
                    
                    contable_sol.lg( false );
                    contable_res.lg( false );
                    contable_num.lg( false );

                    }

                
                
            } 
        
        }
        
        
        
        
        LOG << "Finished Unit Test: " << TestName << nl;
        
        return 0;
}
