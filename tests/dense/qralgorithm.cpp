
#include "../../basic.hpp"
#include "../../dense/factorization.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/factorization.hpp"
#include "../../dense/factorization.hpp"
// #include "../../dense/scalarfunctions.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Tests for QR Algorithm" << nl;
    
    {
        LOG << "1. QR Iteration" << nl;
    
        const int dim = 4;
        DenseMatrix A(dim,dim);

        A.zeromatrix();
        for( int s = 0; s < dim; s++ )
        for( int t = 0; t < dim; t++ )
            A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
            
        int repetitions = 4000;

        FloatVector D = QRIteration( A );
        
        LOG << "Diagonals:" << D << nl;
        
    }
    
    
    
    {
        LOG << "2. Solving Least-Squares problem with QR factorization" << nl;
    
        const int cols =  4;
        const int rows = 40;
        
        DenseMatrix A(rows,cols);
        A.randommatrix();
        
        FloatVector b( A.getdimout() );
        b.random();
        
        auto x = SolveOverconstrained( A, b );
        
        auto r = b - A * x;
        
        LOG << "residual: " << r.norm() << nl;
        LOG << "orthogonality: " << (Transpose(A) * r).norm() << nl;
        
    }
    
    
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
    
    
}
