
#include "../../basic.hpp"
#include "../../utility/stl.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../vtk/vtkwriter.hpp"


using namespace std;

#include "vtk.testsnippet.cxx"

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for VTK output of Simplicial Mesh (3D)" << nl;
    
    {
        
        const MeshSimplicial3D Mx = FicheraCorner3D();  std::string meshname = "Fichera Corner 3D";
        
        internal_print( Mx, meshname );
        
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 3; c++ ) {
            
                M.uniformrefinement();

                M.shake_interior_vertices();
                
                internal_print( M, meshname );
            
            }
            
        }    
        
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 3; c++ ) {
            
                std::vector<int> refinementedges;
                
                for( int k = 0; k < 3 + M.count_edges() / 10; k++ )
                    refinementedges.push_back( random_integer() % M.count_edges() );
                
                sort_and_remove_duplicates( refinementedges );
                
                M.longest_edge_bisection_recursive( refinementedges );

                internal_print( M, meshname );
            
            }
            
        }    

    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
