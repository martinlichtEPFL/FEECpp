
#include <ostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../utility/convergencetable.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;


#include "vtk.testsnippet.cxx"






int main( int argc, char *argv[] )
{
    LOG << "Output of a few important meshes" << nl;
    
    
    {
        
        const int L = 5;
    
        MeshSimplicial2D M = SphericalSurface2D(L);
            
        LOG << L << ":\t" << M.getShapemeasure() << nl;
        
        internal_print( M, "spherical surface 2D", "sphere" );
        
    }
    
    
    
    {
        
        const int L = 5;
    
        MeshSimplicial2D M = LShapedDomain2D();
            
        LOG << L << ":\t" << M.getShapemeasure() << nl;
        
        internal_print( M, "L shaped domain 2D", "lshaped" );
        
    }
    
    
    
    {
        LOG << "Unit Test for VTK output of Simplicial Mesh" << nl;
        
        const int K = 4;
        const int L = 12;
        
        MeshSimplicial2D M = Halo(K,L);
            
        internal_print( M, "halo 2D", "halo" );

    }
    
    
    
    {
        const int Lmin = 4;
        const int Lmax = 4;
        
        for( int L = Lmin; L <= Lmax; L++ )
        {
            
            MeshSimplicial2D M = UnitDisk(L);
            
            LOG << L << ":\t" << M.getShapemeasure() << nl;
            
            int last_original_vertex = M.count_vertices()-1;
            for( int t = 0; t < 15; t++ ) {
                int target_edge1         = M.get_vertex_firstparent_edge(last_original_vertex);
                int target_edge2         = M.get_vertex_nextparent_edge(last_original_vertex,target_edge1);
                int target_edge3         = M.get_vertex_nextparent_edge(last_original_vertex,target_edge2);
                assert( target_edge1 != Mesh::nullindex && target_edge2 != Mesh::nullindex && target_edge3 != Mesh::nullindex );
                // M.newest_vertex_bisection_recursive( target_edge1 );
                // M.newest_vertex_bisection_recursive( target_edge2 );
                // M.newest_vertex_bisection_recursive( target_edge3 );
                // M.newest_vertex_bisection_recursive( target_edge2 );
                // M.newest_vertex_bisection_recursive( target_edge1 );
                // M.newest_vertex_bisection_recursive( target_edge3 );
                // M.newest_vertex_bisection_recursive( target_edge2 );
                M.newest_vertex_bisection_recursive( random_integer() % M.count_edges() );
            }    

            {
                fstream fs( string("./rounddisk.svg"), std::fstream::out );
                // fs << M.outputTikZ();
                fs << M.outputSVG();
                fs.close();
            }
            
            internal_print( M, "round disk 2D", "rounddisk" );

        }
        
    }
    
    
    
    {
        
        const int Lmin = 3;
        const int Lmax = 9;
        
        MeshSimplicial2D M = Annulus( Lmin, Lmax );
        
        internal_print( M, "annulus 2D", "annulus" );
        
    }
        
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
