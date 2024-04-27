
#include <ostream>
#include <fstream>

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../vtk/vtkwriter.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

#include "vtk.testsnippet.cxx"

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for VTK output of Simplicial Mesh (1D)" << nl;
    
    {
        
        MeshSimplicial1D Mx = StandardInterval1D(); string meshname = string("One-dimensional Test Mesh: ") + getbasename(__FILE__);
        
        internal_print( Mx, meshname );
        
        {
            
            auto M = Mx;
            
            for( int c = 0; c < 6; c++ ) {
            
                M.uniformrefinement();

                // M.shake_interior_vertices(); // only if inner and outer dimension match
                
                internal_print( M, meshname );
            
            }
            
        }    
        
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
