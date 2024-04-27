

/**/

#include "../../basic.hpp"
#include "../../operators/composedoperators.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.flags.hpp"
#include "../../utility/convergencetable.hpp"

#include "../../fem/global.trace.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
        
    LOG << "Unit Test: (2D) degree elevation of interpolation preserves mass" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial2D M = UnitTriangle2D();

    M.automatic_dirichlet_flags();
    
    M.check();

    assert( M.getouterdimension() == 2 );
    
    
    const int r_min = 1;
    
    const int r_max = 3;
    
    const int l_min = 0;
    
    const int l_max = 2;

    const int n = M.getinnerdimension();
    
    const int number_of_samples = 1;
        
    
    Float errors[ n+1 ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];

    for( int l = 0; l < l_min; l++ )
        M.uniformrefinement();

    for( int l = l_min; l <= l_max; l++ ){
        
        for( int k = 0;     k <= n; k++ ) 
        for( int r = r_min; r <= r_max;                 r++ ) 
        {
            LOG << nl;

            LOG << "Level:       " << l_min << " <= " << l << " <= " << l_max                 << nl;
            LOG << "Polydegree:  " << r_min << " <= " << r << " <= " << r_max                 << nl;
            LOG << "Form degree: " <<     0 << " <= " << k << " <= " << M.getinnerdimension() << nl;

            LOG << "# T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;


            LOG << "assemble matrices..." << nl;
            
            const auto inclusion  = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), k, r );
            
            const auto trace      = FEECBrokenTraceMatrix( M, M.getinnerdimension(), k, r, true );

            const auto massmatrix = ( M.getinnerdimension() == k ) ? SparseMatrix(ScalingOperator(0,0.)) : FEECBrokenMassMatrix( M, M.getinnerdimension()-1, k, r );
            
            
            
            // if( k != 0 ) continue;
            // auto field = trace.createinputvector();
            // field.setentries(1.0);
            // auto traces_of_field = trace * field;
            
            // // auto e0 = massmatrix.createinputvector(); e0.setentries(0.);
            // // auto e1 = e0; 
            // // auto e2 = e0; 
            
            // // e0[0] = 1.; 
            // // e1[1] = 1.; 
            // // e2[2] = 1.; 
            // // assert( e0.getdimension() == 3 );
            
            // // LOG << M.get_edge_length(0) << tab << M.getMeasure(1,0) << nl; 
            // // LOG << M.get_edge_length(1) << tab << M.getMeasure(1,1) << nl; 
            // // LOG << M.get_edge_length(2) << tab << M.getMeasure(1,2) << nl; 
            
            // // LOG << e0.norm_sq(massmatrix) << nl;
            // // LOG << e1.norm_sq(massmatrix) << nl;
            // // LOG << e2.norm_sq(massmatrix) << nl;

            // // return 0;

            // auto all_edges = massmatrix.createinputvector();
            // all_edges.setentries(1.0);
            // LOG << "Mass of all edges:    " <<       all_edges.norm(massmatrix)    << nl;
            // LOG << "Mass of all edges sq: " <<       all_edges.norm_sq(massmatrix) << nl;
            // LOG << "Trace of constant sq: " << traces_of_field.norm_sq(massmatrix) << nl;
            // continue;






            
            
            

            LOG << massmatrix.getdimout() << 'x' << massmatrix.getdimin() << space 
                << trace.getdimout() << 'x' << trace.getdimin() << space 
                << inclusion.getdimout() << 'x' << inclusion.getdimin() 
                << nl;
                
            
            errors[k][ l-l_min ][ r-r_min ] = 0.;
            

            

            for( int i = 0; i < number_of_samples; i++ ){

                auto field = inclusion.createinputvector();
                field.random();
                field.normalize();

                auto included_field = inclusion * field;

                auto traces_of_field = trace * included_field;

                assert( field.isfinite()           );
                assert( included_field.isfinite()  );
                assert( traces_of_field.isfinite() );

                if( k != M.getinnerdimension() )
                for( int e = 0; e < M.count_simplices(M.getinnerdimension()-1); e++ )
                {
                    const int dim = SullivanSpanSize( n-1, k, r );

                    assert( dim * M.count_simplices(n-1) == traces_of_field.getdimension() );

                    if( M.get_flag( n-1, e ) == SimplexFlag::SimplexFlagDirichlet ) {
                        
                        for( int j = e * dim; j < (e+1)*dim; j++ )
                            assert( traces_of_field[j] == 0. );
                        
                    } else {
                        
                        for( int j = e * dim; j < (e+1)*dim; j++ )
                            assert( is_numerically_small( traces_of_field[j] ) );
                        
                        // LOG << "Non-Dirichlet edge " << e << " with values:" << tab;
                        // for( int j = e * dim; j < (e+1)*dim; j++ )
                        //     LOG << traces_of_field[j] << tab << tab;
                        // LOG << nl;
                        
                    }

                }
                
                const auto error_eucl = traces_of_field.norm();
                
                const auto error_mass = traces_of_field.norm(massmatrix);

                Assert( error_eucl >= -desired_closeness , error_eucl );
                Assert( error_mass >= -desired_closeness , error_mass );

                LOG << traces_of_field.getdimension() << " dimensional with mass " << error_mass << nl; 
                
                Float error = error_mass;

                error = maximum( error, error_eucl );
                
                errors[k][l-l_min][r-r_min] = maximum( errors[k][l-l_min][r-r_min], error );
                
            }
            
        }
        
        if( l != l_max )
        {
            LOG << "Refinement..." << nl;
        
            M.uniformrefinement();

            LOG << "Distortion..." << nl;

            M.shake_interior_vertices();
        }
        
    } 

    LOG << "Convergence tables" << nl;


    std::vector<ConvergenceTable> contables( n+1 );
    
    for( int k = 0; k <= n; k++ ) 
        contables[k].table_name = "Rounding errors D2K" + std::to_string(k);
    for( int k = 0; k <= n; k++ ) {
        for( int r = r_min; r <= r_max; r++ ) 
            contables[k] << printf_into_string("R%d", r-r_min );;
        contables[k] << nl;
    }

    for( int k = 0; k <= n; k++ ) 
    for( int l = l_min; l <= l_max; l++ ) 
    {
        
        for( int r = r_min; r <= r_max; r++ ) 
            contables[k] << errors[k][l-l_min][r-r_min];
        
        contables[k] << nl; 
        
    }
    
    LOG << "Check that differences are small: " << desired_closeness << nl;
    
    for( int k = 0; k <= n; k++ ) 
    {

        contables[k].lg(); 
        LOG << "-------------------" << nl;

        for( int l = l_min; l <= l_max; l++ ) 
        for( int r = r_min; r <= r_max; r++ ) 
        {
            Assert( errors[k][l-l_min][r-r_min] < desired_closeness, errors[k][l-l_min][r-r_min], desired_closeness );
        }

    }
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
