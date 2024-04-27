

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples1D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/lagrangematrices.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (?D) Lagrange matrices agree with FEEC analogues" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial1D M1 = UnitInterval1D();
        MeshSimplicial2D M2 = UnitTriangle2D(); //StandardSquare2D_simple();
        MeshSimplicial3D M3 = UnitSimplex3D();
        
        M1.check();
        M2.check();
        M3.check();
        
        Mesh* Ms[3] = { &M1, &M2, &M3 };
        
        
        
        const int number_of_samples = 50;
        
        const int number_of_comparisons = 9;
        
        const int l_min = 0;
        
        const int l_max = 3;
        
        Float errors[l_max-l_min+1][3][number_of_comparisons];
        
        
        
        for( int l = 0; l < l_min; l++ ) {
            M1.uniformrefinement();
            M2.uniformrefinement();
            M3.uniformrefinement();            
        }
        
        
        
        for( int l = l_min; l <= l_max; l++ )
        {
        
            for( int d = 0; d < 3; d++ )
            {
                
                int m = l - l_min;
                
                for( int t = 0; t < number_of_comparisons; t++ )
                    errors[m][d][t] = -10.;
                
                Mesh& M = *(Ms[d]);
                
                
                LOG << "DIMENSION " << d+1 << " AT LEVEL " << l << nl;
        
                LOG << "...basic FEEC matrices" << nl;
        
                auto feec_broken_mass = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, 1 );

                auto feec_vectormass = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, 0 );
                
                auto feec_diff = FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, 1 );

                auto feec_diff_t = feec_diff.getTranspose();

                auto feec_inc = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, 1 );
                
                auto feec_inc_t = feec_inc.getTranspose();
            
                auto feec_whitney_inc = FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, 1 );
                
                auto feec_whitney_inc_t = feec_whitney_inc.getTranspose();
            
                assert( feec_broken_mass.isfinite() );
                assert( feec_vectormass.isfinite() );
                assert( feec_diff.isfinite() );
                assert( feec_diff_t.isfinite() );
                assert( feec_inc.isfinite() );
                assert( feec_inc_t.isfinite() );
                assert( feec_whitney_inc.isfinite() );
                assert( feec_whitney_inc_t.isfinite() );
                    
                LOG << "...composed FEEC matrices" << nl;
                
                auto feec_broken_stiffness = feec_diff_t & feec_vectormass       & feec_diff;
                auto feec_stiffness        = feec_inc_t  & feec_broken_stiffness & feec_inc;
                auto feec_mass             = feec_inc_t  & feec_broken_mass      & feec_inc;
                
                auto feec_whitney_stiffness        = feec_whitney_inc_t  & feec_broken_stiffness & feec_whitney_inc;
                auto feec_whitney_mass             = feec_whitney_inc_t  & feec_broken_mass      & feec_whitney_inc;
                
                assert( feec_broken_stiffness.isfinite() );
                
                assert( feec_stiffness.isfinite()        );
                assert( feec_mass.isfinite()             );
                
                assert( feec_whitney_stiffness.isfinite() );
                assert( feec_whitney_mass.isfinite()      );
                
                LOG << "...basic Lagrange matrices" << nl;
                
                auto lagr_broken_mass      = LagrangeBrokenMassMatrix( M, 1 );
                auto lagr_broken_stiffness = LagrangeBrokenStiffnessMatrix( M, 1 );
                
                auto lagr_mass      = LagrangeMassMatrix( M, 1 );
                auto lagr_stiffness = LagrangeStiffnessMatrix( M, 1 );
                
                auto lagr_inc   = LagrangeInclusionMatrix( M, M.getinnerdimension(), 1 );
                auto lagr_inc_t = lagr_inc.getTranspose();
                
                LOG << "...composed Lagrange matrices" << nl;
                
                auto lagr_composed_mass      = lagr_inc_t & lagr_broken_mass      & lagr_inc;
                auto lagr_composed_stiffness = lagr_inc_t & lagr_broken_stiffness & lagr_inc;
                
                LOG << "...COMPARISONS" << nl;
                    
                for( int n = 0; n < number_of_samples; n++ ){
                    auto vec = lagr_inc.createinputvector();
                    vec.zero();
                    vec.random();
                    vec.normalize();
                    assert( vec.isfinite() );
                    
                    // inclusion matrices
                    {
                        auto vec_error = ( ( lagr_inc - feec_inc ) * vec ).norm();
                    
                        errors[m][d][0] = maximum( vec_error, errors[m][d][0] );
                    }
                    
                    /*mass matrices (Sullivan) */
                    {
                        auto vec_error = ( ( lagr_mass - feec_mass ) * vec ).norm();
                    
                        errors[m][d][1] = maximum( vec_error, errors[m][d][1] );
                    }
                    
                    /*mass matrices (Whitney) */
                    {
                        auto vec_error = ( ( lagr_mass - feec_whitney_mass ) * vec ).norm();
                    
                        errors[m][d][2] = maximum( vec_error, errors[m][d][1] );
                    }
                    
                    /*mass matrices*/
                    {
                        auto vec_error = ( ( lagr_mass - lagr_composed_mass ) * vec ).norm();
                    
                        errors[m][d][3] = maximum( vec_error, errors[m][d][2] );
                    }
                    
                    /*stiffness matrices (Sullivan)*/
                    {
                        auto vec_error = ( ( lagr_composed_stiffness - feec_stiffness ) * vec ).norm();
                    
                        errors[m][d][4] = maximum( vec_error, errors[m][d][3] );
                    }
                    
                    /*stiffness matrices (Whitney) */
                    {
                        auto vec_error = ( ( lagr_composed_stiffness - feec_whitney_stiffness ) * vec ).norm();
                    
                        errors[m][d][5] = maximum( vec_error, errors[m][d][3] );
                    }
                    
                    /*stiffness matrices */
                    {
                        auto vec_error = ( ( lagr_stiffness - lagr_composed_stiffness ) * vec ).norm();
                    
                        errors[m][d][6] = maximum( vec_error, errors[m][d][4] );
                    }
                    
                }
                
                for( int n = 0; n < number_of_samples; n++ ){
                    auto vec = feec_broken_mass.createinputvector();
                    vec.zero();
                    vec.random();
                    vec.normalize();
                    assert( vec.isfinite() );
                    
                    /*broken mass*/
                    {
                        auto vec_error = ( ( lagr_broken_mass - feec_broken_mass ) * vec ).norm();
                    
                        errors[m][d][7] = maximum( vec_error, errors[m][d][5] );
                    }
                    
                    /*broken stiffness */
                    {
                        auto vec_error = ( ( lagr_broken_stiffness - feec_broken_stiffness ) * vec ).norm();
                    
                        errors[m][d][8] = maximum( vec_error, errors[m][d][6] );
                    }
                    
                }
                
                
                
//                 if( d==2 and l==0 ){
//                     lagr_inc.lg();
//                     lagr_stiffness.lg();
//                     lagr_broken_stiffness.lg();
//                     feec_stiffness.lg();
//                     exit(0);
//                 }
                
                
                
            } // loop over d  
    
                
            if( l != l_max ) { 

                LOG << "Refinement..." << nl;
        
                M1.uniformrefinement();
                M2.uniformrefinement();
                M3.uniformrefinement();
        
                M1.shake_interior_vertices();
                M2.shake_interior_vertices();
                M3.shake_interior_vertices();

            }
        
        }
        
        
        
        {
            
            LOG << "Convergence tables, final results" << nl;

            ConvergenceTable contables[3];
            
            
            
            for( int d = 0; d <            3; d++ )
            {
                
                contables[d].table_name = "Rounding errors, D" + std::to_string(d+1);
                contables[d] << "inc";           // 0
                contables[d] << "mass S";        // 1
                contables[d] << "mass W";        // 2 
                contables[d] << "mass comp";     // 3
                contables[d] << "stiff S";       // 4
                contables[d] << "stiff W";       // 5
                contables[d] << "stiff comp";    // 6
                contables[d] << "br mass";       // 7
                contables[d] << "br stiff";      // 8
                contables[d] << nl;
                
                
                for( int m = 0; m <= l_max-l_min; m++ ) 
                {
                    
                    for( int t = 0; t < number_of_comparisons; t++ )
                    {
                        contables[d] << errors[m][d][t];
                    }
                    
                    contables[d] << nl; 
                    
                }    
            }

            for( int d = 0; d < 3; d++ ) {
                LOG << "Dimension: " << d+1 << '\n';
                contables[d].lg();
                LOG << "----------------------------------" << nl;
            }

        }
            
            
        LOG << "Check that differences are small: " << desired_closeness << nl;
        
        for( int l = l_min; l <= l_max; l++ ) 
        for( int d = 0; d < 3; d++ )
        for( int t = 0; t < number_of_comparisons; t++ )
        {
            if( not ( errors[l-l_min][d][t] < desired_closeness ) )
                LOG << l << space << d << space << t << space << errors[l-l_min][d][t] << nl;
            Assert( errors[l-l_min][d][t] < desired_closeness, errors[l-l_min][d][t], desired_closeness );
        }
        
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
