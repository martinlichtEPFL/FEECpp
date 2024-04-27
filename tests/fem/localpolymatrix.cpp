

/**/

#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"
#include "../../dense/factorization.hpp"
#include "../../mesh/coordinates.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
    LOG << "Unit Test for Inverse of Poly Matrix" << nl;
    
    // LOG << std::setprecision(10);

    const int r_min = 1;
    const int r_max = 9;

    const int n_min = 1;
    const int n_max = 3;

    for( int n = n_min; n <= n_max; n++ )
    for( int r = r_min; r <  r_max; r++ )
    {
        
        DenseMatrix MM = polynomialmassmatrix( n, r );
    
        const int N = MM.getdimin();

        LOG << "Dimension: " << space << n_min << " <= " << n << " <= " << n_max << nl;
        LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

                
        LOG << "Matrix dimension: " << N << nl;
        LOG << "Determinant: " << Determinant(MM) << nl;

        LOG << "Inverse..." << nl;
        DenseMatrix MMinv = Inverse(MM);

        LOG << "Cholesky decomposition..." << nl;
        DenseMatrix MMchol = CholeskyDecomposition(MM);
        
        LOG << "QR decomposition..." << nl;
        DenseMatrix MMqr_q(MM), MMqr_r(MM);
        QRFactorization(MM,MMqr_q,MMqr_r);
        
        const Float diff_inv  = ( MM * MMinv - IdentityMatrix(N) ).norm();
        const Float diff_chol = ( MMchol * Transpose(MMchol) - MM ).norm();
        const Float diff_qr   = ( MMqr_q * MMqr_r - MM ).norm();

        LOG << "\ta=" << diff_inv
            << "\tb=" << diff_chol
            << "\tc=" << diff_qr
            << nl;
                    
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
    
    return 0;
}
