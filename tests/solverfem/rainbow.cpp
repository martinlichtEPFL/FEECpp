
#include "../../basic.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../sparse/matcsr.hpp"
#include "../../sparse/rainbow.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.massmatrix.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit test: (3D) CSR-matrix rainbow partition" << nl;
    
    LOG << "Initial mesh..." << nl;
    
    MeshSimplicial3D M = StandardCube3D();
    
    M.check();

    
    
    
    const int r = 1;
    
    const int l = 4;
    
    for( int i = 0; i < l; i++ ) M.uniformrefinement();

    
    LOG << "Assemble matrices..." << nl;
    
    SparseMatrix incmatrix_scalar   = FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r );
    SparseMatrix incmatrix_scalar_t = incmatrix_scalar.getTranspose();
    SparseMatrix massmatrix_scalar  = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );

    MatrixCSR mass = MatrixCSR( incmatrix_scalar_t & massmatrix_scalar & incmatrix_scalar );

    LOG << "Rainbow..." << nl;
        
    const Rainbow rainbow( mass );

    LOG << "Number of rows:    " << mass.getdimout() << nl;
    LOG << "Maximum row width: " << mass.getmaxrowwidth() << nl;
    LOG << "Number of colors:  " << rainbow.num_colors << nl;
    LOG << "Rough Lower bound: " << binomial_integer( 3+r, r ) << nl;

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
