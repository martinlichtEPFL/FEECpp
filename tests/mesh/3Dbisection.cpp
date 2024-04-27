
#include "../../basic.hpp"
#include "../../utility/random.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 3D Module" << nl;
    
    for( int ei = 0; ei < 6; ei++ )
    {
        
        MeshSimplicial3D M = UnitSimplex3D();
        
        M.check();
        
        M.automatic_dirichlet_flags();

        M.check_dirichlet_flags();
        
        for( int c = 0; c < 10; c++ ) {
            M.bisect_edge( M.get_tetrahedron_edge( 0, ei ) );
        }
        
        M.check();
        
        M.check_dirichlet_flags();
        
    }
    
    
    {
        
        MeshSimplicial3D M = UnitCube3D();
        
        M.check();
        
        M.automatic_dirichlet_flags();

        for( int c = 0; c < 20; c++ ) {
            M.bisect_edge( M.get_tetrahedron_edge( c % 5, random_integer() % 6 ) );
        }
        
        M.check();
        
        M.check_dirichlet_flags();

    }
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
