

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.coefficientmassmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test for Weighted Mass Matrices" << nl;
    
    LOG << "Initial mesh..." << nl;
        
    MeshSimplicial3D M = SomeSimplex3D();
    for( int i = 0; i < 8; i++ )
        M.bisect_edge(i);
    
    
    M.check();

    std::function<DenseMatrix(const FloatVector&)> generator
            = [](const FloatVector&) -> DenseMatrix{ return IdentityMatrix(3); };
    
    const int min_r = 0;
    const int max_r = 4;
    
    const int min_s = 0;
    const int max_s = 4;
    
    for( int r = min_r; r <= max_r; r++ ) 
    for( int s = min_s; s <= max_s; s++ ) 
    {
        
        LOG << "r=" << r << " s=" << s << nl;

        LOG << "vector mass matrix, classical" << nl;

        SparseMatrix mass1 = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
        
        LOG << "vector mass matrix, weighted" << nl;

        SparseMatrix mass1_s = FEECBrokenCoefficientMassMatrix( 
                                        M, M.getinnerdimension(), 1, r, s, generator );
        
        LOG << "pseudo mass matrix, classical" << nl;

        SparseMatrix mass2 = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
        
        LOG << "pseudo mass matrix, weighted" << nl;

        SparseMatrix mass2_s = FEECBrokenCoefficientMassMatrix( 
                                        M, M.getinnerdimension(), 2, r, s, generator );
        
        LOG << "comparing products with random vectors" << nl;
        
        if( s==1 and r==0)
            LOG << mass1 << nl << mass1_s << nl << mass2 << nl << mass2_s << nl;

        const int L = 6;
        for( int i = 0; i < L; i++ )
        {
            FloatVector v1 = mass1.createinputvector();
            FloatVector v2 = mass2.createinputvector();
            v1.random(); v1.normalize(mass1);
            v2.random(); v2.normalize(mass2);

            auto diff1 = mass1 * v1 - mass1_s * v1;
            auto diff2 = mass2 * v2 - mass2_s * v2;

            Float norm_diff1 = diff1.norm();
            Float norm_diff2 = diff2.norm();

            Assert( norm_diff1 < desired_closeness, norm_diff1, desired_closeness );
            Assert( norm_diff2 < desired_closeness, norm_diff2, desired_closeness );
        }
        LOG << nl;
        
    }


    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
