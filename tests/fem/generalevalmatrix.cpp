

/**/

#include "../../basic.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/utilities.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        LOG << "Unit Test: Evaluation Matrix and its Invertibility" << nl;
        
        // LOG << std::setprecision(10);
        
        const int n_min = 1;
        const int n_max = 3;
        const int r_min = 1;
        const int r_max = 6;

        for( int n = n_min; n <= n_max; n++ )
        for( int r = r_min; r <= r_max; r++ )
        {
            LOG << "Dimension: " << space << n_min << " <= " << n << " <= " << n_max << nl;
            LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;
            
            const auto lpsbc = InterpolationPointsInBarycentricCoordinates( n, r );
            
            const auto EM = PointValuesOfMonomials( r, lpsbc );
            
            assert( EM.issquare() );
        
            const auto EMinv = Inverse( EM );
        
            int N = EM.getdimin();

            Float diff_inv_1 = ( EM * EMinv - IdentityMatrix(N) ).norm();
            Float diff_inv_2 = ( EMinv * EM - IdentityMatrix(N) ).norm();
            
            LOG << "dim(A)=" << N << "\tdiff1=" << diff_inv_1 << "\tdiff2=" << diff_inv_2 << space << machine_epsilon << nl;
            
            assert( diff_inv_1 < desired_closeness );
            assert( diff_inv_2 < desired_closeness );
            
        }
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        
        return 0;
}
