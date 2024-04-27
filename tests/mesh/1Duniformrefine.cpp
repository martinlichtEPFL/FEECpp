
#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial1D.hpp"
#include "../../mesh/examples1D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Manifold 1D Module" << nl;

    MeshSimplicial1D M = StandardInterval1D();
    
    LOG << "Check" << nl;
    
    M.check();

    LOG << M << nl;
    
    LOG << "Start refinement" << nl;
    
    for( int c = 0; c <= 10; c++ )
    {
        
        LOG << "uniform refinement: " << c << nl;

        M.improved_uniformrefinement();
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
