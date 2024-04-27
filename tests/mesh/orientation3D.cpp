
#include <vector>

#include "../../basic.hpp"
#include "../../combinatorics/generateindexmaps.hpp"
#include "../../utility/stl.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Simplicial 3D Module" << nl;
    
    
    MeshSimplicial3D M = UnitSimplex3D();

    const int l_min = 0;
    const int l_max = 2;

    for( int c = 0; c < l_min; c++ )
        M.uniformrefinement();
    

    auto face_inclusions = generateSigmas( IndexRange(0,2), IndexRange(0,3) );

    for( int c = 0; c <= l_max; c++ )
    {
        
        for( int f = 0; f < M.count_faces(); f++ )
        {
            
            auto ftp = M.get_tetrahedron_parents_of_face(f);

            auto fv = M.getsubsimplices(2,0,f);

            for( auto t : ftp ) 
            {
                
                auto tv = M.getsubsimplices(3,0,t);

                int index = 0;
                for( ; index < face_inclusions.size(); index++ )
                    if( tv * face_inclusions[index] == fv )
                        break;
                assert( index < face_inclusions.size() );

                int gap_index = 0;
                while( gap_index <= M.getinnerdimension()-1 && fv[gap_index] == tv[gap_index] ) gap_index++;

                int geometric_orientation = sign_integer( Determinant( M.getTransformationJacobian(3,t) ) );

                LOGPRINTF( 
                    "Face %3d (%3d %3d %3d) \thas parent %3d (%3d %3d %3d %3d) \twith orientation s * (-1)^%3d = %2d \n", 
                    f, fv[0], fv[1], fv[2],
                    t, tv[0], tv[1], tv[2], tv[3],
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
