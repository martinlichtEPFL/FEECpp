
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"
#include "../../sparse/sparsematrix.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Test for SparseMatrix" << nl;

    {
        LOG << "1. Action of Identity Matrix" << nl;
        const int dim = 5;
        
        SparseMatrix M( dim, dim );

        for( int i = 0; i < dim; i++ ) M.appendentry( i, i, 1. );
        
        FloatVector vec( dim, 1.23 );
        for( int i = 0; i < dim; i++ ) vec[i] = i * 1.2345;

        auto Mvec = M * vec;

        for( int i = 0; i < dim; i++ ) assert( vec[i] == Mvec[i] );
        
    }
    
    {
        SparseMatrix M( 2, 3 );

        for( int i = 0; i < 5; i++ )
            for( int j = 0; j < 7; j++ )
                M.appendentry( (3*i) % 2, (2*j) % 3, i / 3. + j*j );

        LOG << "This is the content of some matrix:" << nl;
        LOG << M << nl;

        LOG << "Sort the Entries" << nl;
        M.sortentries();
        LOG << M << nl;

        M.clearentries();
        LOG << "Empty Matrix again" << nl;
        LOG << M << nl;

        LOG << "Next Matrix:" << nl;
        M.appendentry( 0, 0, 1. );
        M.appendentry( 0, 1, 2. );
        M.appendentry( 0, 2, 3. );
        M.appendentry( 1, 0, 4. );
        M.appendentry( 1, 1, 5. );
        M.appendentry( 1, 2, 6. );
        LOG << M << nl;

        FloatVector vec(3);
        vec.setentry(0,13);
        vec.setentry(1,17);
        vec.setentry(2,19);
        LOG << "Some vector:" << nl;
        LOG << vec << nl;

        LOG << "Matrix-Vector Product:" << nl;
        LOG << M * vec << nl;
    
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
}
