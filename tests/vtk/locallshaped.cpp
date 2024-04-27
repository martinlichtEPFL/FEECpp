
#include <algorithm>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for VTK output of Simplicial Mesh" << nl;
    
    // MeshSimplicial2D M = UnitCubeTriangulation(3,3);
    MeshSimplicial2D M = LShapedDomain2D();
    
    int l_max = 5;

    for( int l = 0; l < 7; l++ )
    {
        LOG << "Print VTK-type file" << nl;
        LOG << "T/E/V: " << M.count_triangles() << "/" << M.count_edges() << "/" << M.count_vertices() << nl;
        
        fstream fs( string("./locallshaped") + std::to_string(l) + string(".vtk"), std::fstream::out );

        VTKWriter vtk( M, fs, "L-Shaped Domain" );
        // vtk.writeCoordinateBlock();
        // vtk.writeTopDimensionalCells();

        fs.close();

        LOG << "Refine" << nl;

        if( l != l_max ) {
            
            FloatVector cellwisemass = FloatVector( M.count_triangles(), 
                [&M]( int t ) -> Float{
                    FloatVector mp = M.get_triangle_midpoint(t);
                    return exp( - sqrt( mp[0]*mp[0] + mp[1]*mp[1] ) );
                }
            );
            
            Float maxcellwisemass = cellwisemass.maxnorm();

            std::vector<int> marked_edges;
            marked_edges.reserve( 3 * M.count_edges() );

            for( int s = 0; s < M.count_triangles(); s++ ) 
            if( cellwisemass.at(s) > 0.75 * maxcellwisemass )
            {
                // LOG << M.get_triangle_edge( s, 0 ) << space << M.get_triangle_edge( s, 1 ) << space << M.get_triangle_edge( s, 2 );\\ << nl;
                marked_edges.push_back( M.get_triangle_edge( s, 0 ) );
                marked_edges.push_back( M.get_triangle_edge( s, 1 ) );
                marked_edges.push_back( M.get_triangle_edge( s, 2 ) );
            }

            std::sort( marked_edges.begin(), marked_edges.end() );
            auto temp = std::unique( marked_edges.begin(), marked_edges.end() );
            marked_edges.erase( temp, marked_edges.end() );
            
            LOG << "marked edges: " << marked_edges.size() << "/" << M.count_edges() << nl;

            // for( int e = 0; e < M.count_edges(); e++ )
            //     if( e % 10 == 0)
            //         marked_edges.push_back( e );
            
            M.newest_vertex_bisection_recursive( marked_edges );

        }


    }
        
        
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
