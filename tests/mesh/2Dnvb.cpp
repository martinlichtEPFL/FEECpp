
#include <vector>

#include "../../basic.hpp"
#include "../../utility/random.hpp"
#include "../../utility/stl.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 2D Module" << nl;
    
    {
        
        LOG << "1. Refine triangles" << nl;
        
        MeshSimplicial2D M = UnitTriangle2D();
        
        M.check();

        M.automatic_dirichlet_flags();

        M.check_dirichlet_flags();
        
        for( int k = 0; k <= 2; k++ ) 
        { 
            LOG << "Uniform refinements..." << k << nl;
            M.uniformrefinement();
        }
        
        LOG << "Newest vertex bisections..." << nl;

        int cell_count_initial = M.count_triangles();
        int cell_marked_count  = 0;
        
        int c_max = 50;
        
        for( int c = 0; c < c_max; c++ ) {
        
            std::vector<int> markedcells;
            
            unsigned int p = 40;
            for( int t = 0; t < M.count_triangles(); t++ )
                if( random_integer() % p == 0 ) 
                    markedcells.push_back( t );
            markedcells.push_back( 0 );
            cell_marked_count += markedcells.size();
            
            std::vector<int> markededges;
            for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
            sort_and_remove_duplicates( markededges );
            
            LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();

        M.check_dirichlet_flags();
        
    }
    
    
    {
        
        LOG << "2. Refine surface, choose cells via uniform distribution" << nl;
        
        MeshSimplicial2D M = TetrahedralSurface2D(); 
        
        M.check();

        M.automatic_dirichlet_flags();

        M.check_dirichlet_flags();
        
        for( int k = 0; k <= 2; k++ ) 
        { 
            LOG << "Uniform refinements..." << k << nl;
            M.uniformrefinement();
        }
        
        LOG << "Newest vertex bisections..." << nl;

        int cell_count_initial = M.count_triangles();
        int cell_marked_count  = 0;
        
        int c_max = 20;
        
        for( int c = 0; c < c_max; c++ ) {
        
            std::vector<int> markedcells;
            
            unsigned int p = 10;
            for( int t = 0; t < M.count_triangles(); t++ )
                if( random_integer() % p == 0 ) 
                    markedcells.push_back( t );
            markedcells.push_back( 0 );
            cell_marked_count += markedcells.size();
            
            std::vector<int> markededges;
            for( int t : markedcells ) markededges.push_back( M.get_oldest_edge( t ) );
            sort_and_remove_duplicates( markededges );
            
            LOG << c << "/" << c_max << " Refine " << markedcells.size() << "/" << M.count_triangles() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_triangles() - cell_count_initial )/(Float)( cell_marked_count ) << nl;
        
        }
        
        M.check();

        M.check_dirichlet_flags();
        
    }
    
    {
        
        LOG << "3. NVB, choose edges via uniform distribution" << nl;
        
        MeshSimplicial2D M = UnitTriangle2D();
        
        M.check();

        M.automatic_dirichlet_flags();

        M.check_dirichlet_flags();
        
        LOG << "Newest vertex bisections..." << nl;

        int edges_count_initial = M.count_edges();
        int edges_marked_count  = 0;
        
        int iter_max = 4000;
        
        for( int i = 0; i < iter_max; i++ ) {
        
            std::vector<int> markededges;
            markededges.push_back( random_integer() % M.count_edges() );
            sort_and_remove_duplicates( markededges );

            edges_marked_count += markededges.size();
            
            LOG << i << "/" << iter_max << " Refine " << markededges.size() << "/" << M.count_edges() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_edges() - edges_count_initial )/(Float)( edges_marked_count ) << nl;
        
        }
        
        M.check();

        M.check_dirichlet_flags();
        
    }
    
    
    {
        
        LOG << "4. repeated bisection of a fixed triangle" << nl;
        
        MeshSimplicial2D M = UnitTriangle2D();
        
        M.check();

        M.automatic_dirichlet_flags();

        M.check_dirichlet_flags();
        
        int edges_count_initial = M.count_edges();
        int edges_marked_count  = 0;
        
        int iter_max = 2000;
        
        LOG << "Newest vertex bisections..." << nl;

        for( int i = 0; i < iter_max; i++ ) {
        
            std::vector<int> markededges;
            markededges.push_back( 0 );
            sort_and_remove_duplicates( markededges );

            edges_marked_count += markededges.size();
            
            LOG << i << "/" << iter_max << " Refine " << markededges.size() << "/" << M.count_edges() << " ... ";
            M.newest_vertex_bisection_recursive( markededges );
            LOG << "Ratio=" << ( M.count_edges() - edges_count_initial )/(Float)( edges_marked_count ) << nl;
        
        }
        
        M.check();

        M.check_dirichlet_flags();
        
    }


    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
