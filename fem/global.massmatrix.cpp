
#include <vector>

#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/factorization.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/polynomialmassmatrix.hpp"

#include "../fem/global.massmatrix.hpp"



inline DenseMatrix calcAtA( const DenseMatrix& A ) { return Transpose(A) * A; }

SparseMatrix FEECBrokenMassMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    if(false) 
    {
        DenseMatrix multiplier( n, n+1, 0. );
        for( int i = 0; i <  n; i++ ) {
            multiplier(i,i+1) =  1.;
            multiplier(i,  0) = -1.;
        }
        
        DenseMatrix rankone = IdentityMatrix(n+1) - DenseMatrix( n+1,n+1, 1./(n+1) );
        
        LOG << multiplier << nl;
        LOG << rankone << nl;
        LOG << multiplier * rankone << nl;

    }
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
    
    const int localdim = binomial_integer( n + r, n ) * binomial_integer( n+1, k );
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r );

    assert( polyMM.issquare() );
    assert( polyMM.getdimin() == binomial_integer( n+r, n ) );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        Float measure      = mesh.getMeasure( n, s );

        assert( measure >= 0. );

        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
        // DenseMatrix GPM    = calcAtA( mesh.getGradientMatrix( n, s ) );
            
        assert( GPM.isfinite() );
        if( not ( GPM - calcAtA( mesh.getGradientMatrix( n, s ) ) ).is_numerically_small() ){
            const auto AtA = calcAtA( mesh.getGradientMatrix( n, s ) );
            assert( AtA.isfinite() );
            LOG << AtA << nl;
            LOG << GPM << nl;
            LOG << ( GPM - AtA ).norm() << nl;
            assert( ( GPM - AtA ).is_numerically_small() );
        }

        if(false)
        {
            auto result = GPM * FloatVector( n+1, 1. );

            LOG << result.norm() << nl;

            Assert( result.is_numerically_small(), result );
        }

        if(false)
        {
            DenseMatrix rankone = IdentityMatrix(n+1) - DenseMatrix( n+1,n+1, 1./(n+1) );

            const DenseMatrix foo = Transpose(rankone) * GPM * rankone;

            const DenseMatrix delta = GPM - foo;

            Assert( delta.is_numerically_small(), GPM, foo );
            
            // LOGPRINTF("%.15Le\n", (long double)delta.sumnorm() );

            // GPM = foo;
        }

        DenseMatrix formMM = SubdeterminantMatrix( GPM, k );


        if(false)
        {
            const DenseMatrix Jac = mesh.getTransformationJacobian( n, s );

            DenseMatrix multiplier( n, n+1, 0. );
            for( int i = 0; i <  n; i++ ) {
                multiplier(i,i+1) =  1.;
                multiplier(i,  0) = -1.;
            }

            if(false) // TODO: this leads to failure ...
            {
                const DenseMatrix middle = Inverse( Transpose(Jac) * Jac );

                const DenseMatrix middle_minors = SubdeterminantMatrix( middle, k );

                const DenseMatrix multiplier_minors = SubdeterminantMatrix( multiplier, k );

                // LOG << multiplier_minors << nl;

                const DenseMatrix foo = Transpose(multiplier_minors) * middle_minors * multiplier_minors;

                const DenseMatrix delta = formMM - foo;

                Assert( delta.is_numerically_small(), formMM, foo );

                {
                    auto D = QRIteration( formMM );
                    LOG << k << ":NULL:"<< D << nl;
                }
                {
                    auto D = QRIteration( foo );
                    LOG << k << ":EINS:"<< D << nl;
                }
                // LOGPRINTF("%.15Le\n", (long double)delta.sumnorm() );

                // formMM = foo;
            }
            
            if(false)
            {
                DenseMatrix R( Jac.getdimin() );
                DenseMatrix Q( Jac.getdimout(), Jac.getdimin() );
                QRFactorization( Jac, Q, R );

                const DenseMatrix Rinv = Inverse(R);
                const DenseMatrix Rinvt = Transpose(Rinv);
        
                const DenseMatrix minors = SubdeterminantMatrix( Rinvt * multiplier, k );

                const DenseMatrix foo = Transpose(minors) * minors;

                const DenseMatrix delta = formMM - foo;

                Assert( delta.is_numerically_small(), formMM, foo );

                {
                    auto D = QRIteration( foo );
                    LOG << k << ":ZWEI:"<< D << nl;
                }
                // LOGPRINTF("%.15Le\n", (long double)delta.sumnorm() );
                
                // formMM = foo;
            }

            if(false)
            {
                DenseMatrix rankone = IdentityMatrix(n+1) - DenseMatrix( n+1,n+1, 1./(n+1) );

                DenseMatrix rankone_minors = SubdeterminantMatrix(rankone,k);

                const DenseMatrix foo = Transpose(rankone_minors) * formMM * rankone_minors;

                const DenseMatrix delta = formMM - foo;

                Assert( delta.is_numerically_small(), formMM, foo );

                {
                    auto D = QRIteration( foo );
                    LOG << k << ":DREI:"<< D << nl; // INDICATES WORSE PERFORMANCE
                }
                // LOGPRINTF("%.15Le\n", (long double)delta.sumnorm() );

                // formMM = foo;
            }
        }
        
        
        if(false)
        { // This does have an effect somewhere ...
            DenseMatrix Aux1( n+1, n+1, 0. );
            for( int i = 1; i <= n; i++ ) {
                Aux1(i,i) = 1.;
                Aux1(i,0) = -1.;
            }
            
            const DenseMatrix Aux2 = SubdeterminantMatrix( Aux1, k );

            const DenseMatrix foo = formMM;

            {
                auto D = QRIteration( foo );
                LOG << k << ":VIER:" << D << nl;
            }
            // formMM = Transpose(Aux2) * formMM * Aux2;

            const DenseMatrix delta = formMM - foo;

            Assert( delta.is_numerically_small(), formMM, foo );
        }

        // {
        //     auto result = formMM * FloatVector( formMM.getdimout(), 1. );
        //     LOG << result.norm() << nl;
        //     Assert( result.is_numerically_small(), result );
        // }

        // LOG << formMM << nl;

        
        
        DenseMatrix fullMM = MatrixTensorProduct( polyMM, formMM ) * measure;

        if(false)
        {
            DenseMatrix Aux1 = IdentityMatrix(n+1) - DenseMatrix( n+1, n+1, 1./(n+1) );
            const DenseMatrix Aux2 = SubdeterminantMatrix( Aux1, k );
            const DenseMatrix Aux3 = DenseMatrix( polyMM.getdimin(), Aux2, 1. );
            auto foo = Transpose(Aux3) * fullMM * Aux3; 
            fullMM = foo;
        }
        
        
        if(false)
        { // This does have an effect anywhere ...
            DenseMatrix Aux1( n+1, n+1, 0. );
            for( int i = 1; i <= n; i++ ) {
                Aux1(i,i) = 1.;
                Aux1(i,0) = -1.;
            }
            
            const DenseMatrix Aux2 = MatrixTensorProduct( IdentityMatrix(polyMM.getdimin()), SubdeterminantMatrix( Aux1, k ) );

            const DenseMatrix foo = fullMM;

            // fullMM = Transpose(Aux2) * fullMM * Aux2;

            const DenseMatrix delta = fullMM - foo;

            Assert( delta.is_numerically_small(), fullMM, foo );

            {
                auto D = QRIteration( fullMM );
                LOG << k << ":FUNF:"<< D << nl;
            }
            {
                auto D = QRIteration( foo );
                LOG << k << ":SECHS:"<< D << nl;
            }
            // LOGPRINTF( "delta %.10e\n", delta.maxnorm() );
        }
        
        assert( k != 0 or ( fullMM - polyMM * measure ).is_numerically_small() );
        
        for( int i = 0; i < localdim; i++ )
        for( int j = 0; j < localdim; j++ )
        {
            int index_of_entry = s * localdim * localdim + i * localdim + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + i;
            entry.column = s * localdim + j;
            entry.value  = fullMM( i, j );
            
            ret.setentry( index_of_entry, entry );
        }
        
    }
    
    return ret;
}





SparseMatrix FEECBrokenMassMatrixRightFactor( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim_in  = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    const int localdim_out = binomial_integer( n+r, n ) * binomial_integer( n  , k );

    const int dim_in      = num_simplices * localdim_in;
    const int dim_out     = num_simplices * localdim_out;
    const int num_entries = num_simplices * localdim_in * localdim_out;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r );
    
    DenseMatrix polyMM_right = Transpose(CholeskyDecomposition(polyMM));

    {

        // LOG << polyMM << nl;        
        // LOG << ( Transpose(polyMM_right) * polyMM_right ) << nl;        
        // LOG << polyMM_right.getdimin() << nl;        
        // LOG << polyMM.getdimin() << space << ( Transpose(polyMM_right) * polyMM_right - polyMM ).norm() << nl;        
        assert( ( Transpose(polyMM_right) * polyMM_right - polyMM ).is_numerically_small() ); 

    }
            
//     LOG << polyMM << nl;
        
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        Float measure = mesh.getMeasure( n, s );
            
        DenseMatrix GPM_right    = mesh.getGradientProductMatrixRightFactor( n, s );
            
        DenseMatrix formMM_right = SubdeterminantMatrix( GPM_right, k );
    
        DenseMatrix fullMM_right = MatrixTensorProduct( polyMM_right, formMM_right ) * std::sqrt(measure);

        assert( fullMM_right.getdimin()  == formMM_right.getdimin()  * polyMM_right.getdimin()  );
        assert( fullMM_right.getdimout() == formMM_right.getdimout() * polyMM_right.getdimout() );
        
        {
            
            /* CHECK WHETHER THE ARITHMETICS WORK OUT */
            
            DenseMatrix polyMM_0 = polynomialmassmatrix( n, r );
            DenseMatrix formMM_0 = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, s ), k );
            DenseMatrix fullMM_0 = MatrixTensorProduct( polyMM_0, formMM_0 ) * measure;
            
            assert( ( Transpose(formMM_right) * formMM_right - formMM_0 ).is_numerically_small() ); 
            assert( ( Transpose(fullMM_right) * fullMM_right - fullMM_0 ).is_numerically_small() ); 
            
        }
        
        
        for( int i = 0; i < localdim_out; i++ )
        for( int j = 0; j < localdim_in ; j++ )
        {
            int index_of_entry = s * localdim_out * localdim_in + i * localdim_in + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim_out + i;
            entry.column = s * localdim_in  + j;
            entry.value  = fullMM_right( i, j );
            
            ret.setentry( index_of_entry, entry );
        }
        
        
        
    }
    
    return ret;
}





FloatVector FEECBrokenMassMatrix_cellwisemass( const Mesh& mesh, int n, int k, int r, const FloatVector vec )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
        
    FloatVector ret( num_simplices );
    
    assert( vec.getdimension() == localdim * num_simplices );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        ret[s] = 0.;
        for( int i = 0; i < localdim; i++ )
            ret[s] = ret[s] + vec[ s * localdim + i] * vec[ s * localdim + i];
        
    }
    
    return ret;
}



SparseMatrix FEECBrokenMassMatrix_cellwiseinverse( const Mesh& mesh, int n, int k, int r )
{
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    DenseMatrix polyMM = polynomialmassmatrix( n, r );
    
    DenseMatrix polyMMinv = Inverse( polyMM ); 
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        Float measure      = mesh.getMeasure( n, s );

        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
        
        auto aux1   = GPM + DenseMatrix( n+1, n+1, 1. );
        auto aux2   = Inverse( aux1 );
        auto GPMinv = aux2 - DenseMatrix( n+1, n+1, 1. );
            
        DenseMatrix formMMinv = SubdeterminantMatrix( GPMinv, k );
    
        DenseMatrix fullMMinv = MatrixTensorProduct( polyMMinv, formMMinv ) / measure;

        for( int i = 0; i < localdim; i++ )
        for( int j = 0; j < localdim; j++ )
        {
            int index_of_entry = s * localdim * localdim + i * localdim + j;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + i;
            entry.column = s * localdim + j;
            entry.value  = fullMMinv( i, j );
            
            ret.setentry( index_of_entry, entry );
        }
        
    }
    
    return ret;
}








// inline DenseMatrix elementmassmatrix( int n, int ambientdim, int r, int k, DenseMatrix Jacobian )
// {
//     
//     assert( Jacobian.getdimin() == n && Jacobian.getdimout() == ambientdim ); 
//     assert( n <= ambientdim );
//     assert( n >= 0 && n >= k );
//         
//     DenseMatrix polyMM = polynomialmassmatrix( n, r );
//     
//     DenseMatrix formMM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, t ), k );
//         
//     return TensorProduct( polyMM, formMM ) * absolute( determinant( Jacobian ) ) / factorial_numerical( n );
//     
// }

