
#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 2D Module" << nl;
    
    MeshSimplicial2D M = StandardSquare2D();
    
    M.check();

    M.automatic_dirichlet_flags();

    M.check_dirichlet_flags();
    
    LOG << "Refinement" << nl;
    
    for( int c = 0; c < 5; c++ ) {

        M.midpoint_refinement_global();
        
    }
        
    M.check();

    M.check_dirichlet_flags();
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
