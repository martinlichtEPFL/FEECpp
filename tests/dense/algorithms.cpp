
#include "../../basic.hpp"
#include "../../dense/factorization.hpp"
#include "../../dense/functions.hpp"
#include "../../dense/factorization.hpp"
#include "../../dense/simplesolver.hpp"
// #include "../../dense/scalarfunctions.hpp"


int main( int argc, char *argv[] )
{
    LOG << "Unit Tests for Matrix Algorithms" << nl;
    
    // LOG << std::setprecision(5);
    // LOG << std::fixed << std::ios::floatfield; //setf( std::ios::fixed, std::ios::floatfield ); 
    // LOG << std::showpos;
    
    {
        LOG << "A. Scalar functions of matrices" << nl;
       
        DenseMatrix S( 4, 4 );
        S(0,0) =  3; S(0,1) =  0; S(0,2) = 6; S(0,3) =  0; 
        S(1,0) =  1; S(1,1) = -1; S(1,2) = 0; S(1,3) =  2; 
        S(2,0) = -1; S(2,1) =  1; S(2,2) = 1; S(2,3) = -1; 
        S(3,0) =  2; S(3,1) = -4; S(3,2) = 4; S(3,3) =  0; 
        
        
        LOG << S << nl;
        LOG << "Matrix trace:   " << S.trace() << nl;
        LOG << "Norm L1:        " << S.sumnorm() << nl;
        LOG << "Norm Frobenius: " << S.frobeniusnorm() << nl;
        LOG << "Norm Max:       " << S.maxnorm() << nl;
        
        Float p = 1.01;
        LOG << nl << "Norm Lp with p=" << p << ": " << S.lpnorm( p ) << nl << nl;
        
        Float p1 = 100.0001;
        Float p2 = 1.00001;
        LOG << nl << "Row " << p1 << space
                     << "Col " << p2 << space
                     << S.norm_row_col( p1, p2 ) << nl;
        LOG << nl << "Col " << p1 << space
                     << "Row " << p2 << space 
                     << S.norm_col_row( p1, p2 ) << nl;
        
        LOG << nl << "Row " << 1. << space
                     << "Col " << 1. << space
                     << S.norm_row_col( 1., 1. ) << nl;
        LOG << nl << "Col " << 1. << space
                     << "Row " << 1. << space 
                     << S.norm_col_row( 1., 1. ) << nl;
        LOG << nl;
        
        LOG << nl << "Row " << 2. << space
                     << "Col " << 2. << space
                     << S.norm_row_col( 2., 2. ) << nl;
        LOG << nl << "Col " << 2. << space
                     << "Row " << 2. << space 
                     << S.norm_col_row( 2., 2. ) << nl;
        LOG << nl;
        
        LOG << nl << "Row " << 20. << space
                     << "Col " << 20. << space
                     << S.norm_row_col( 20., 20. ) << nl;
        LOG << nl << "Col " << 20. << space
                     << "Row " << 20. << space 
                     << S.norm_col_row( 20., 20. ) << nl;
        LOG << nl;
        
        
        
        LOG << "Norm Operator L1:  " << S.NormOperatorL1() << nl;
        LOG << "Norm Operator Max: " << S.NormOperatorMax() << nl;
        LOG << nl;
        
        LOG << "GerschgorinRow:    " << S.GerschgorinRow() << nl;
        LOG << "GerschgorinColumn: " << S.GerschgorinColumn() << nl;
        
    }
    
    
    {
        LOG << "B. Transposing dense matrices" << nl;
        
        LOG << "Transpose of square matrices" << nl;
        
        for( int dim = 0; dim < 4; dim++ ) 
        {
            
            std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                        { return 2*r + c; };
            
            DenseMatrix A(dim, testmatrix);
            DenseMatrix B = A;
            
            auto At1 = Transpose( A );
            auto At2 = TransposeSquare( A );
            LOG << A << At1 << At2 << nl;

            TransposeInSitu( A );
            TransposeSquareInSitu( B );
            LOG << A << B << nl;
            
        }
        
        LOG << "Transpose of non-square matrices" << nl;
        
        for( int dimr = 0; dimr < 4; dimr++ ) 
        for( int dimc = 0; dimc < 4; dimc++ ) 
        {
            
            std::function<Float(int,int)> testmatrix = [](int r, int c) 
                                                        { return 2*r + c; };
            
            DenseMatrix A( dimr, dimc, testmatrix );
            
            LOG << A << Transpose(A) << nl;
            TransposeInSitu( A );
            LOG << A << nl;
            
        }
        
    }
    
    
    {
        
        LOG << "1. Cholesky decomposition" << nl;
        
        DenseMatrix A(3,3);
        
        A(0,0) =   4; A(0,1) =  12; A(0,2) = -16; 
        A(1,0) =  12; A(1,1) =  37; A(1,2) = -43; 
        A(2,0) = -16; A(2,1) = -43; A(2,2) =  98; 
        
        DenseMatrix L = CholeskyDecomposition( A );

        DenseMatrix ActualL( 3, 3 );
        ActualL(0,0) =  2; ActualL(0,1) = 0; ActualL(0,2) = 0; 
        ActualL(1,0) =  6; ActualL(1,1) = 1; ActualL(1,2) = 0; 
        ActualL(2,0) = -8; ActualL(2,1) = 5; ActualL(2,2) = 3; 
        
        //LOG << "Original matrix:" << A << nl;
        //LOG << "Factor matrix:" << L << nl;
        //LOG << "Product matrix:" << L * Transpose(L) << nl;

        assert( ( ActualL - L ).is_numerically_small() );
        assert( ( L * Transpose(L) - A ).is_numerically_small() );
        
    }
      
    {
    
        LOG << "2. Cholesky decomposition" << nl;
        
        int dim = 3;
        DenseMatrix A(dim);
        
        A.zeromatrix();
        for( int s = 0; s < dim; s++ )
        for( int t = 0; t < dim; t++ )
            A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
        
        DenseMatrix L = CholeskyDecomposition( A );
        
        LOG << "Original matrix:" << A << nl;
        
        LOG << "Factor matrix:" << L << nl;
        
        LOG << "Product matrix:" << L * Transpose(L) << nl;

        assert( ( L * Transpose(L) - A ).is_numerically_small() );
        
    }
    
    {
        LOG << "3. Matrix determinant, inverse, cofactor matrix" << nl;
        
        for( int dim = 1; dim < 6; dim++ ) 
        {
            DenseMatrix A = HilbertMatrix( dim );
            
            LOG << A << nl;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );

            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;
            LOG << 1. / Determinant(A) - Determinant( Ainv ) << nl;
            LOG << A - Inverse( Ainv ) << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );
        }
    }
    
    {
      
        LOG << "4. Matrix determinant, inverse, cofactor matrix" << nl;

        DenseMatrix A( 3 );
        A(0,0) = -2; A(0,1) =  2; A(0,2) = -3; 
        A(1,0) = -1; A(1,1) =  1; A(1,2) =  3; 
        A(2,0) =  2; A(2,1) =  0; A(2,2) = -1; 

        Float det = Determinant(A);
        Float det_l = Determinant_laplaceexpansion(A);
        Float det_g = Determinant_gauss(A);
        Float det_b = Determinant_bareiss(A);

        LOG << "Determinant (default): " << det << nl;
        LOG << "Determinant (laplace): " << det_l << nl;
        LOG << "Determinant (gauss):   " << det_g << nl;
        LOG << "Determinant (bareiss): " << det_b << nl;

        assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );

        const auto Acof = CofactorMatrix( A );
        const auto Ainv = Inverse( A );

        LOG << Acof << nl;
        LOG << Ainv << nl;
        LOG << A * Ainv << nl;
        LOG << Ainv * A << nl;

        assert( ( Ainv * A ).is_numerically_identity() );
        assert( ( A * Ainv ).is_numerically_identity() );
      
    }
    
    {
        
        LOG << "5. Gauss Jordan algorithm" << nl;
    
        DenseMatrix A(4,4);
        
        A(0,0) = 1; A(0,1) = 2; A(0,2) = 3; A(0,3) = 5; 
        A(1,0) = 1; A(1,1) = 1; A(1,2) = 1; A(1,3) = 6; 
        A(2,0) = 3; A(2,1) = 3; A(2,2) = 1; A(2,3) = 2; 
        A(3,0) = 2; A(3,1) = 1; A(3,2) = 0; A(3,3) = 1; 
        
        LOG << "Determinant (default): " << Determinant(A) << nl;
        LOG << "Determinant (laplace): " << Determinant_laplaceexpansion(A) << nl;
        LOG << "Determinant (gauss):   " << Determinant_gauss(A) << nl;
      
        
//         A(0,0) = 4; A(0,1) = 1; A(0,2) = 0; 
//         A(1,0) = 0; A(1,1) = 1; A(1,2) = 0; 
//         A(2,0) = 1; A(2,1) = 0; A(2,2) = 1; 
        
        LOG << "Original matrix:" << A << nl;
        
        const DenseMatrix Ainv = GaussJordanInplace( A );
        
        LOG << "Proposed Inverse:" << Ainv << nl;
        
        LOG << "Product of both:" << Ainv * A << nl;

        assert( ( Ainv * A ).is_numerically_identity() );
        
    }
    
    {
        LOG << "6. Gauss Jordan algorithm" << nl;
        for( int i = 0; i < 6; i++ )
        {
            
            DenseMatrix A(8);
            LOG << A << nl;
            A.randomintegermatrix(-2,2);
            LOG << A << nl;
            
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;
            
            LOG << A << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );

            const auto Ainv = GaussJordanInplace(A);
            
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
        }
    }

    {
        LOG << "7. Gauss Jordan algorithm" << nl;
    
        for( int N = 1; N <= 7; N++ )
        {
            DenseMatrix C(N);
            for( int i = 0; i < N; i++ )
            for( int j = 0; j < N; j++ )
                C(i,j) = 1. / ( i+j+1 );
            
            LOG << C * GaussJordanInplace(C,true) << nl;
            
            {
                auto Cinv = GaussJordanInplace(C);
                const auto I = IdentityMatrix(N);
                // Float relaxation = Cinv.norm();

                Assert( ( C * Cinv ).is_numerically_identity(), C * Cinv );
            }
        }
        
        
    }
    
    
    {
        LOG << "8. QR and QL Factorization" << nl;
    
        for( int dim = 1; dim <= 10; dim++ ){
            
            DenseMatrix A(dim,dim);
        
            A.zeromatrix();
            for( int s = 0; s < dim; s++ )
            for( int t = 0; t < dim; t++ )
                A(s,t) = 3 * kronecker(s,t) - kronecker(s,t-1) - kronecker(s,t+1);
                
            {
                DenseMatrix Q(dim,dim), R(dim,dim);
                
                QRFactorization( A, Q, R );
                
                DenseMatrix Rinv =   Inverse(R);
                DenseMatrix Qt   = Transpose(Q);
                
                assert( ( Q * R - A ).is_numerically_small() );
                assert( ( Q * Qt ).is_numerically_identity() );
                
                LOG << "Matrix A:" << A;
                LOG << "Matrix Q:" << Q;
                LOG << "Matrix R:" << R;
                LOG << "Matrix Q * R:" << Q * R;
                LOG << "Matrix Q^t:" << Qt;
                LOG << "Matrix Rinv:" << Rinv;
                LOG << "Matrix R * Rinv:" << R * Rinv;
                LOG << "Matrix Rinv * R:" << Rinv * R;
                LOG << "Matrix Q * Qinv:" << Q * Qt;
                LOG << "Matrix Qinv * Q:" << Qt * Q;
                LOG << "Matrix Rinv * Q^t * A:" << Rinv * Qt * A;
                LOG << "Matrix A * Rinv * Q^t:" << A * Rinv * Qt;
            }

            {
                DenseMatrix L(dim,dim), Q(dim,dim);
                
                LQFactorization( A, L, Q );
                
                DenseMatrix Linv =   Inverse(L);
                DenseMatrix Qt   = Transpose(Q);

                assert( ( L * Q - A ).is_numerically_small() );
                assert( ( Q * Qt ).is_numerically_identity() );
                
                LOG << "Matrix A:" << A;
                LOG << "Matrix Q:" << Q;
                LOG << "Matrix L:" << L;
                LOG << "Matrix L * Q:" << L * Q;
                LOG << "Matrix Q^t:" << Qt;
                LOG << "Matrix Linv:" << Linv;
                LOG << "Matrix L * Linv:" << L * Linv;
                LOG << "Matrix Linv * L:" << Linv * L;
                LOG << "Matrix Q * Qinv:" << Q * Qt;
                LOG << "Matrix Qinv * Q:" << Qt * Q;
                LOG << "Matrix Q^t * Linv * A:" << Qt * Linv * A;
                LOG << "Matrix A * Qt * Linv:" << A * Qt * Linv;
            }

            {
                DenseMatrix Q1(dim,dim), R(dim,dim);
                QRFactorization( A, Q1, R );

                DenseMatrix L(dim,dim), Q2(dim,dim);
                LQFactorization( Transpose(A), L, Q2 );

                assert( ( Q1 * R - Transpose(Q2) * Transpose(L) ).is_numerically_small() );
                
            }


        }
        
        
    }
    
    // if(false)
    {
        LOG << "9. Compare Determinats of Matrices with random coefficients" << nl;
        for( int t = 0; t < 7; t++ )
        for( int i = 0; i < 6; i++ )
        {
            
            DenseMatrix A(t);
            A.randomintegermatrix(10,80);
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;

            // LOGPRINTF( "%.17e %.17e %.17e %.17e \n", det, det_l, det_g, det_b );

            assert( is_numerically_close( det / det_l, 1. ) && is_numerically_close( det / det_g, 1. ) && is_numerically_close( det / det_b, 1. ) );
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );
        }
    }
    
    {
        LOG << "10. Inverse and determinants of unit matrices" << nl;
        for( int t = 0; t < 7; t++ )
        {
            
            DenseMatrix A = 3. * IdentityMatrix(t);
            
            // LOG << A << nl;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = GaussJordanInplace( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );
        }
    }

    
    // if(false)
    {
        LOG << "11. Compare Determinats of random orthogonal Matrices" << nl;
        for( int t = 3; t < 7; t++ )
        for( int i = 0; i < 6; i++ )
        {
            
            DenseMatrix A(t);
            A.random_orthogonal_matrix();

            LOG << A << nl;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << "Determinant (default): " << det << nl;
            LOG << "Determinant (laplace): " << det_l << nl;
            LOG << "Determinant (gauss):   " << det_g << nl;
            LOG << "Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );
            
            {
                const int dim = t;
                DenseMatrix Q(dim,dim), R(dim,dim);
                QRFactorization( A, Q, R );
                Float detR = UpperTriangularDeterminant(R); Float detQ = Determinant(Q);
                LOG << "Determinant (QR): " << detR << nl;
                Assert( is_numerically_close( detQ, 1. ) or is_numerically_close( detQ, -1. ));
                Assert( is_numerically_close( det, detQ * detR ) );

            }
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );
        }
    }
    
    
    {
        LOG << "12. Compare Determinats of random ill-conditioned 2x2 matrices" << nl;
        for( int t = 3; t < 5; t++ )
        {
            
            DenseMatrix A(2,2);
            A(0,0) = power_numerical( 10, t );
            A(0,1) = power_numerical( 10, t ) - 1;
            A(1,0) = power_numerical( 10, t ) - 1;
            A(1,1) = power_numerical( 10, t ) - 2;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << t << " Determinant (default): " << det << nl;
            LOG << t << " Determinant (laplace): " << det_l << nl;
            LOG << t << " Determinant (gauss):   " << det_g << nl;
            LOG << t << " Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );
            
            {
                const int dim = 2;
                DenseMatrix Q(dim,dim), R(dim,dim);
                QRFactorization( A, Q, R );
                Float detR = UpperTriangularDeterminant(R); Float detQ = Determinant(Q);
                LOG << "Determinant (QR): " << detR << nl;
                Assert( is_numerically_close( detQ, 1. ) or is_numerically_close( detQ, -1. ));
                Assert( is_numerically_close( det, detQ * detR ) );
            }
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity( 1e-7 ), Ainv * A );
            assert( ( A * Ainv ).is_numerically_identity( 1e-7 ), A * Ainv );
        }
    }
    

    
    
    {
        LOG << "13. Compare Determinats of random ill-conditioned matrices" << nl;
        for( int t = 3; t < 5; t++ )
        {
            
            DenseMatrix Q(2);
            Q.random_orthogonal_matrix();

            DenseMatrix M(2,2);
            M(0,0) = power_numerical( 10, t );
            M(0,1) = 0.;
            M(1,0) = 0.;
            M(1,1) = 1. / power_numerical( 10, t );

            DenseMatrix A = Transpose(Q) * M * Q;
            
            Float det = Determinant(A);
            Float det_l = Determinant_laplaceexpansion(A);
            Float det_g = Determinant_gauss(A);
            Float det_b = Determinant_bareiss(A);

            LOG << t << " Determinant (default): " << det << nl;
            LOG << t << " Determinant (laplace): " << det_l << nl;
            LOG << t << " Determinant (gauss):   " << det_g << nl;
            LOG << t << " Determinant (bareiss): " << det_b << nl;

            assert( is_numerically_close( det, det_l ) && is_numerically_close( det, det_g ) && is_numerically_close( det, det_b ) );
            
            {
                const int dim = 2;
                DenseMatrix Q(dim,dim), R(dim,dim);
                QRFactorization( A, Q, R );
                Float detR = UpperTriangularDeterminant(R); Float detQ = Determinant(Q);
                LOG << "Determinant (QR): " << detR << nl;
                Assert( is_numerically_close( detQ, 1. ) or is_numerically_close( detQ, -1. ));
                Assert( is_numerically_close( det, detQ * detR ) );
            }
            
            const auto Acof = CofactorMatrix( A );
            const auto Ainv = Inverse( A );

            LOG << Acof << nl;
            LOG << Ainv << nl;
            LOG << A * Ainv << nl;
            LOG << Ainv * A << nl;

            assert( ( Ainv * A ).is_numerically_identity() );
            assert( ( A * Ainv ).is_numerically_identity() );
        }
    }
    
    LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

    return 0;
    
    
}
