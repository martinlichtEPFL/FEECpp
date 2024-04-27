
#include <sstream>

#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/io.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 3D Module" << nl;

    MeshSimplicial3D M = UnitSimplex3D();

    M.check();

    LOG << "Set automatic Dirichlet flags..." << nl;

    M.automatic_dirichlet_flags();

    M.check();

    M.check_dirichlet_flags();

    LOG << "Refinement..." << nl;

    for( int c = 0; c < 3; c++ ) 
    {
        M.uniformrefinement();

        LOG << c << "\tpatch size:          " << M.getVertexPatchSize() << nl;
        LOG << c << "\tmaximum diameter:    " << M.getMaximumDiameter() << nl;
        LOG << c << "\tminumum diameter:    " << M.getMinimumDiameter() << nl;
        LOG << c << "\tcomparison quotient: " << M.getComparisonQuotient() << nl;
        LOG << c << "\tradii quotient:      " << M.getRadiiQuotient(3) << nl;
        LOG << c << "\theight ratio:        " << M.getHeightQuotient() << nl;
        LOG << c << "\tshape measure:       " << M.getShapemeasure() << nl;
        LOG << nl;
    }

    {
        
        LOG << "start IO..." << nl;
        
        std::stringstream ss;
        
        writeMeshSimplicial3D( ss, M );
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicial3D M2 = readMeshSimplicial3D( ss );
        
        LOG << "check mesh equivalence..." << nl;
        assert( M == M2 );
    }

    M.check();

    M.check_dirichlet_flags();

    {

        for( int e = 0; e < M.count_edges(); e++ )
            assert( is_numerically_close( M.get_edge_length(e), M.getMeasure(1,e) ) );

    }

    LOG << "Standard output..." << nl;

    //LOG << M << nl;

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
