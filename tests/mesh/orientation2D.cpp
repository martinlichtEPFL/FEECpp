
#include <vector>

#include "../../basic.hpp"
#include "../../combinatorics/generateindexmaps.hpp"
#include "../../utility/stl.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 2D Module" << nl;
    
    
    MeshSimplicial2D M = UnitTriangle2D();

    const int l_min = 0;
    const int l_max = 2;

    for( int c = 0; c < l_min; c++ )
        M.uniformrefinement();
    

    auto edge_inclusions = generateSigmas( IndexRange(0,1), IndexRange(0,2) );

    for( int c = 0; c <= 2; c++ )
    {
        
        for( int e = 0; e < M.count_edges(); e++ )
        {
            
            auto etp = M.get_triangle_parents_of_edge(e);

            auto ev = M.getsubsimplices(1,0,e);

            for( auto t : etp ) 
            {
                
                auto tv = M.getsubsimplices(2,0,t);

                int index = 0;
                for( ; index < edge_inclusions.size(); index++ )
                    if( tv * edge_inclusions[index] == ev )
                        break;
                assert( index < edge_inclusions.size() );

                int gap_index = 0;
                while( gap_index <= M.getinnerdimension()-1 && ev[gap_index] == tv[gap_index] ) gap_index++;

                int geometric_orientation = sign_integer( Determinant( M.getTransformationJacobian(2,t) ) );

                LOGPRINTF( "Edge %3d (%3d %3d) \thas parent %3d (%3d %3d %3d) \twith orientation s * (-1)^%3d = %2d \n", 
                    e, ev[0], ev[1],
                    t, tv[0], tv[1], tv[2],
                    gap_index, sign_power( gap_index ) * geometric_orientation
                    );

            }

        }

        LOG << nl;

        if( c != l_max ) M.uniformrefinement();

    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
