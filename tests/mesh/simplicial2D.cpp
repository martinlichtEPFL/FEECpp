
#include <cstdio>

#include <fstream>
#include <sstream>

#include "../../basic.hpp"
#include "../../utility/files.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/io.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 2D Module" << nl;

    // MeshSimplicial2D M = TetrahedralSurface2D();
    MeshSimplicial2D M = StandardSquare2D_strange14();

    M.check();

    LOG << "Set automatic Dirichlet flags..." << nl;

    M.automatic_dirichlet_flags();

    M.check();

    M.check_dirichlet_flags();

    LOG << "Refinement" << nl;

    for( int c = 0; c < 6; c++ ) 
    {
        M.uniformrefinement();

        LOG << c << "\tpatch size:          " << M.getVertexPatchSize() << nl;
        LOG << c << "\tmaximum diameter:    " << M.getMaximumDiameter() << nl;
        LOG << c << "\tminumum diameter:    " << M.getMinimumDiameter() << nl;
        LOG << c << "\tcomparison quotient: " << M.getComparisonQuotient() << nl;
        LOG << c << "\tradii quotient:      " << M.getRadiiQuotient(2) << nl;
        LOG << c << "\theight ratio:        " << M.getHeightQuotient() << nl;
        LOG << c << "\tshape measure:       " << M.getShapemeasure() << nl;
        LOG << nl;
    }

    M.check();

    M.check_dirichlet_flags();

    {
        
        LOG << "start IO..." << nl;
        
        std::stringstream ss;
        
        writeMeshSimplicial2D( ss, M );
        
        ss.seekg( std::ios_base::beg );
        
        MeshSimplicial2D M2 = readMeshSimplicial2D( ss );
        
        LOG << "check mesh equivalence..." << nl;
        assert( M == M2 );
    }

    {

        for( int e = 0; e < M.count_edges(); e++ )
            assert( is_numerically_close( M.get_edge_length(e), M.getMeasure(1,e) ) );

    }

    // LOG << "Standard output..." << nl;
    // LOG << M << nl;

    {
        LOG << "TikZ output..." << nl;

        fstream fs( experimentfile( getbasename(__FILE__), "tex" ), std::fstream::out );
        
        fs << M.outputTikZ();

        fs.close();
    }
    
    {
        LOG << "SVG output..." << nl;

        fstream fs( experimentfile( getbasename(__FILE__), "svg" ), std::fstream::out );
        int num_tets = M.count_triangles();
        FloatVector red( num_tets, 128 ), green( num_tets, 240 ), blue( num_tets, 38 );
        red.random_within_range(0.,255.); green.random_within_range(0.,255.); blue.random_within_range(0.,255.);
        fs << M.outputSVG( 0.01, "array", "white", &red, &green, &blue );
        fs.close();
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
