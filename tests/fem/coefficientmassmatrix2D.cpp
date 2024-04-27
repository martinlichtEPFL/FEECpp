

/**/

#include <fstream>

#include "../../basic.hpp"
#include "../../sparse/sparsematrix.hpp"
#include "../../mesh/mesh.simplicial2D.hpp"
#include "../../mesh/examples2D.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.coefficientmassmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    
    LOG << "Unit Test for Weighted Mass Matrices" << nl;
    
    LOG << "Initial mesh..." << nl;
        
    MeshSimplicial2D M = StandardSquare2D_simple();

    for( int l = 0; l < 5; l++ )
    {
        M.uniformrefinement();
    }
                
    M.shake_interior_vertices();
    
    M.check();

    std::function<DenseMatrix(const FloatVector&)> generator
            = [](const FloatVector&) -> DenseMatrix{ return IdentityMatrix(2); };
    
    const int min_r = 0;
    const int max_r = 4;
    
    const int min_s = 0;
    const int max_s = 4;
    
    for( int r = min_r; r <= max_r; r++ ) 
    for( int s = min_s; s <= max_s; s++ ) 
    {
        
        LOG << "r=" << r << " s=" << s << nl;

        LOG << "vector mass matrix, classical" << nl;

        SparseMatrix massmatrix = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
        
        LOG << "vector mass matrix, weighted" << nl;

        SparseMatrix massmatrix_s = FEECBrokenCoefficientMassMatrix( 
                                        M, M.getinnerdimension(), 1, r, s, generator );
        
        LOG << "comparing products with random vectors" << nl;
        
        assert( massmatrix.getdimin()  == massmatrix_s.getdimin()  );
        assert( massmatrix.getdimout() == massmatrix_s.getdimout() );
        
        const int L = 6;
        for( int i = 0; i < L; i++ )
        {
            FloatVector v = massmatrix.createinputvector();
            v.random();
            v.normalize(massmatrix);

            auto diff = massmatrix * v - massmatrix_s * v;

            Float norm_diff = diff.norm();

            Assert( norm_diff < desired_closeness, norm_diff, desired_closeness );
        }
        LOG << nl;

    }

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
