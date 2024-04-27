
#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Dense Matrix class" << nl;

    LOG << "Random matrix A of size 3x4, and 3*A" << nl;

    DenseMatrix A( 3, 4 );
    A.randommatrix();

    LOG << A << nl;
    LOG << 3 * A << nl;

    LOG << "Random matrix B of size 3x4, and A+B" << nl;

    DenseMatrix B( 3, 4 );
    B.randommatrix();

    LOG << B << nl;
    LOG << A + B << nl;

    LOG << "Unit Matrices of size 3x3 and 4x4" << nl;

    DenseMatrix I3(3,3);
    I3.unitmatrix();
    DenseMatrix I4(4,4);
    I4.unitmatrix();
    LOG << I3 << I4 << nl;
    
    LOG << "I3 * A and A * I4" << nl;

    LOG << I3 * A << nl;
    LOG << A * I4 << nl;
    
    LOG << "5 * I3 and (5*I3)* A" << nl;

    auto S5 = 5. * I3;
    LOG << S5 << nl;
    LOG << operator*( S5, A ) << nl;

    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
