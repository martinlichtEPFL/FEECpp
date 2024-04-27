

/**/

#include "../../basic.hpp"
#include "../../dense/densematrix.hpp"
#include "../../dense/factorization.hpp"
#include "../../dense/factorization.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.diffmatrix.hpp"
#include "../../fem/global.sullivanincl.hpp"
#include "../../fem/global.whitneyincl.hpp"
#include "../../fem/utilities.hpp"
#include "../../solver/iterativesolver.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

// TODO: replace by pivoted row and column elimination

static Float lowest( const FloatVector vec, Float threshold )
{
    assert( vec.getdimension() > 0 and threshold >= -desired_closeness );
    Float ret = vec[0];
    for( int p = 1; p < vec.getdimension(); p++ )
        if( vec[p] < ret and vec[p] > threshold ) ret = vec[p];
    assert( ret >= threshold );
    return ret;
}

static Float PowerMethod( const LinearOperator& A, int repetitions )
{
    assert( A.issquare() );

    auto x = FloatVector( A.getdimin() );
    x.random(); x.normalize();

    while ( repetitions --> 0 )
    {
        x = A * x; x.normalize();
    }
    
    return ( x * (A * x) ) / (x*x);
}


UNUSED static Float InversePowerMethod( const LinearOperator& A, int repetitions )
{
    assert( A.issquare() );

    auto x = FloatVector( A.getdimin() );
    x.random(); x.normalize();

    ConjugateResidualMethod CRM(A);
    CRM.verbosity = IterativeSolver::VerbosityLevel::silent;

    while ( repetitions --> 0 )
    {
        auto y = FloatVector(A.getdimin(),0.);
        CRM.solve(y,x);
        y.normalize();
        x = y;
    }
    
    return ( x * (A * x) ) / (x*x);
}


// static FloatVector Eigenvalues( DenseMatrix A, int repetitions, Float shift = 0. )
// {
//     assert( A.issquare() );
//     const int dim = A.getdimin();
//     DenseMatrix Q(dim,dim), R(dim,dim);
//     for( int p = 0; p < dim; p++) A(p,p) += shift;
//     while ( repetitions --> 0 )
//     {
//         QRFactorization( A, Q, R );
//         A = R * Q;
//     }
//     FloatVector ret(dim);
//     for( int p = 0; p < dim; p++) ret[p] = A(p,p);
//     return ret;        
// }

static FloatVector Eigenvalues( DenseMatrix A, int repetitions, Float shift = 0. )
{
    assert( A.issquare() );
    const int dim = A.getdimin();
    DenseMatrix L(dim,dim);
    for( int p = 0; p < dim; p++) A(p,p) += shift;
    while ( repetitions --> 0 )
    {
        L = CholeskyDecomposition( A );
        A = Transpose(L) * L;
    }
    FloatVector ret(dim);
    for( int p = 0; p < dim; p++) ret[p] = A(p,p);
    return ret;        
}

// static FloatVector Eigenvalues( DenseMatrix A, int repetitions, Float shift = 0. )
// {
//     assert( A.issquare() );
//     const int dim = A.getdimin();
//     for( int p = 0; p < dim; p++) A(p,p) += shift;
    
//     for( int i = 0; i < dim; i++ )
//     for( int r1 =    0; r1 < dim; r1++ )
//     for( int r2 = r1+1; r2 < dim; r2++ )
//     {
//         Float alpha = - A(r2,r1) / A(r1,r1);
//         for( int t = 0; t < dim; t++ ) A(r2,t) = A(r2,t) + alpha * A(r1,t); //A.addrow(c,r,alpha); A.addcolumn(c,r,alpha);
//     }

//     FloatVector ret(dim);
//     for( int p = 0; p < dim; p++) ret[p] = A(p,p);
//     return ret;        
// }

// static FloatVector Eigenvalues( DenseMatrix A, int repetitions, Float shift = 0. )
// {
//     assert( A.issquare() );
//     const int dim = A.getdimin();
//     for( int p = 0; p < dim; p++) A(p,p) += shift;
//     DenseMatrix X(dim,dim), Q(dim,dim), R(dim,dim);
//     X.randommatrix();
//     while ( repetitions --> 0 )
//     {
//         X = A * X;
//         QRFactorization(X,Q,R);
//         X = Q;
//     }
//     FloatVector ret(dim);
//     for( int p = 0; p < dim; p++) ret[p] = R(p,p);
//     return ret;        
// }
    
int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (3D) condition numbers" << nl;
        
        MeshSimplicial3D M = RegularSimplex3D();
        
        M.check();
        
        const int r_min = 1;
        
        const int r_max = 2;
        
        Float min_s_scalar_mass[ r_max - r_min + 1 ];
        Float min_s_vector_mass[ r_max - r_min + 1 ];
        Float min_s_pseudo_mass[ r_max - r_min + 1 ];
        Float min_s_volume_mass[ r_max - r_min + 1 ];
        
        Float max_s_scalar_mass[ r_max - r_min + 1 ];
        Float max_s_vector_mass[ r_max - r_min + 1 ];
        Float max_s_pseudo_mass[ r_max - r_min + 1 ];
        Float max_s_volume_mass[ r_max - r_min + 1 ];
        
        Float min_w_scalar_mass[ r_max - r_min + 1 ];
        Float min_w_vector_mass[ r_max - r_min + 1 ];
        Float min_w_pseudo_mass[ r_max - r_min + 1 ];
        Float min_w_volume_mass[ r_max - r_min + 1 ];
        
        Float max_w_scalar_mass[ r_max - r_min + 1 ];
        Float max_w_vector_mass[ r_max - r_min + 1 ];
        Float max_w_pseudo_mass[ r_max - r_min + 1 ];
        Float max_w_volume_mass[ r_max - r_min + 1 ];
        
        Float min_s_scalar_stiff[ r_max - r_min + 1 ];
        Float min_s_vector_stiff[ r_max - r_min + 1 ];
        Float min_s_pseudo_stiff[ r_max - r_min + 1 ];
        
        Float max_s_scalar_stiff[ r_max - r_min + 1 ];
        Float max_s_vector_stiff[ r_max - r_min + 1 ];
        Float max_s_pseudo_stiff[ r_max - r_min + 1 ];
        
        Float min_w_scalar_stiff[ r_max - r_min + 1 ];
        Float min_w_vector_stiff[ r_max - r_min + 1 ];
        Float min_w_pseudo_stiff[ r_max - r_min + 1 ];
        
        Float max_w_scalar_stiff[ r_max - r_min + 1 ];
        Float max_w_vector_stiff[ r_max - r_min + 1 ];
        Float max_w_pseudo_stiff[ r_max - r_min + 1 ];
        
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;

                LOG << "...assemble mass matrices" << nl;
                
                auto massmatrix_scalar = DenseMatrix( FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r ) );
                auto massmatrix_vector = DenseMatrix( FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r ) );
                auto massmatrix_pseudo = DenseMatrix( FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r ) );  
                auto massmatrix_volume = DenseMatrix( FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r ) );
                
                auto stiffmatrix_scalar = Transpose(DenseMatrix( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r   ) )) *
                                          DenseMatrix( FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r-1 ) ) *
                                          DenseMatrix( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 0, r   ) );

                auto stiffmatrix_vector = Transpose(DenseMatrix( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r   ) )) *
                                          DenseMatrix( FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r-1 ) ) *
                                          DenseMatrix( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 1, r   ) );
                                         
                auto stiffmatrix_pseudo = Transpose(DenseMatrix( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 2, r   ) )) *
                                          DenseMatrix( FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r-1 ) ) *
                                          DenseMatrix( FEECBrokenDiffMatrix( M, M.getinnerdimension(), 2, r   ) );
                                         
                
                auto s_incmatrix_scalar = DenseMatrix( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 0, r ) );
                auto s_incmatrix_vector = DenseMatrix( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 1, r ) );
                auto s_incmatrix_pseudo = DenseMatrix( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 2, r ) );
                auto s_incmatrix_volume = DenseMatrix( FEECSullivanInclusionMatrix( M, M.getinnerdimension(), 3, r ) );

                auto w_incmatrix_scalar = DenseMatrix( FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 0, r ) );
                auto w_incmatrix_vector = DenseMatrix( FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 1, r ) );
                auto w_incmatrix_pseudo = DenseMatrix( FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 2, r ) );
                auto w_incmatrix_volume = DenseMatrix( FEECWhitneyInclusionMatrix( M, M.getinnerdimension(), 3, r ) );

                auto s_mass_scalar = Transpose(s_incmatrix_scalar) * massmatrix_scalar * s_incmatrix_scalar;
                auto s_mass_vector = Transpose(s_incmatrix_vector) * massmatrix_vector * s_incmatrix_vector;
                auto s_mass_pseudo = Transpose(s_incmatrix_pseudo) * massmatrix_pseudo * s_incmatrix_pseudo;
                auto s_mass_volume = Transpose(s_incmatrix_volume) * massmatrix_volume * s_incmatrix_volume;
                
                auto w_mass_scalar = Transpose(w_incmatrix_scalar) * massmatrix_scalar * w_incmatrix_scalar;
                auto w_mass_vector = Transpose(w_incmatrix_vector) * massmatrix_vector * w_incmatrix_vector;
                auto w_mass_pseudo = Transpose(w_incmatrix_pseudo) * massmatrix_pseudo * w_incmatrix_pseudo;
                auto w_mass_volume = Transpose(w_incmatrix_volume) * massmatrix_volume * w_incmatrix_volume;
                
                auto s_stiff_scalar = Transpose(s_incmatrix_scalar) * stiffmatrix_scalar * s_incmatrix_scalar;
                auto s_stiff_vector = Transpose(s_incmatrix_vector) * stiffmatrix_vector * s_incmatrix_vector;
                auto s_stiff_pseudo = Transpose(s_incmatrix_pseudo) * stiffmatrix_pseudo * s_incmatrix_pseudo;
                
                auto w_stiff_scalar = Transpose(w_incmatrix_scalar) * stiffmatrix_scalar * w_incmatrix_scalar;
                auto w_stiff_vector = Transpose(w_incmatrix_vector) * stiffmatrix_vector * w_incmatrix_vector;
                auto w_stiff_pseudo = Transpose(w_incmatrix_pseudo) * stiffmatrix_pseudo * w_incmatrix_pseudo;
                
                assert( s_mass_scalar.isfinite() ); assert( w_mass_scalar.isfinite() ); 
                assert( s_mass_vector.isfinite() ); assert( w_mass_vector.isfinite() );
                assert( s_mass_pseudo.isfinite() ); assert( w_mass_pseudo.isfinite() ); 
                assert( s_mass_volume.isfinite() ); assert( w_mass_volume.isfinite() ); 
                
                assert( s_stiff_scalar.isfinite() ); assert( w_stiff_scalar.isfinite() ); 
                assert( s_stiff_vector.isfinite() ); assert( w_stiff_vector.isfinite() );
                assert( s_stiff_pseudo.isfinite() ); assert( w_stiff_pseudo.isfinite() ); 
                
                LOG << "...diagonalize" << nl;
                
                const int repetitions = 200;
                
                const Float shift = 100 * machine_epsilon;

                auto s_diagonal_mass_scalar = Eigenvalues( s_mass_scalar, repetitions );
                auto s_diagonal_mass_vector = Eigenvalues( s_mass_vector, repetitions );
                auto s_diagonal_mass_pseudo = Eigenvalues( s_mass_pseudo, repetitions );
                auto s_diagonal_mass_volume = Eigenvalues( s_mass_volume, repetitions );
                
                auto w_diagonal_mass_scalar = Eigenvalues( w_mass_scalar, repetitions );
                auto w_diagonal_mass_vector = Eigenvalues( w_mass_vector, repetitions );
                auto w_diagonal_mass_pseudo = Eigenvalues( w_mass_pseudo, repetitions );
                auto w_diagonal_mass_volume = Eigenvalues( w_mass_volume, repetitions );
                
                auto s_diagonal_stiff_scalar = Eigenvalues( s_stiff_scalar, repetitions, shift );
                auto s_diagonal_stiff_vector = Eigenvalues( s_stiff_vector, repetitions, shift );
                auto s_diagonal_stiff_pseudo = Eigenvalues( s_stiff_pseudo, repetitions, shift );
                
                auto w_diagonal_stiff_scalar = Eigenvalues( w_stiff_scalar, repetitions, shift );
                auto w_diagonal_stiff_vector = Eigenvalues( w_stiff_vector, repetitions, shift );
                auto w_diagonal_stiff_pseudo = Eigenvalues( w_stiff_pseudo, repetitions, shift );
                
                assert( s_diagonal_mass_scalar.isfinite() ); assert( s_diagonal_mass_scalar.isnonnegative() );
                assert( s_diagonal_mass_vector.isfinite() ); assert( s_diagonal_mass_vector.isnonnegative() );
                assert( s_diagonal_mass_pseudo.isfinite() ); assert( s_diagonal_mass_pseudo.isnonnegative() );
                assert( s_diagonal_mass_volume.isfinite() ); assert( s_diagonal_mass_volume.isnonnegative() );
                
                assert( w_diagonal_mass_scalar.isfinite() ); assert( w_diagonal_mass_scalar.isnonnegative() );
                assert( w_diagonal_mass_vector.isfinite() ); assert( w_diagonal_mass_vector.isnonnegative() );
                assert( w_diagonal_mass_pseudo.isfinite() ); assert( w_diagonal_mass_pseudo.isnonnegative() );
                assert( w_diagonal_mass_volume.isfinite() ); assert( w_diagonal_mass_volume.isnonnegative() );

                assert( s_diagonal_stiff_scalar.isfinite() ); assert( s_diagonal_stiff_scalar.isnonnegative() );
                assert( s_diagonal_stiff_vector.isfinite() ); assert( s_diagonal_stiff_vector.isnonnegative() );
                assert( s_diagonal_stiff_pseudo.isfinite() ); assert( s_diagonal_stiff_pseudo.isnonnegative() );
                
                assert( w_diagonal_stiff_scalar.isfinite() ); assert( w_diagonal_stiff_scalar.isnonnegative() );
                assert( w_diagonal_stiff_vector.isfinite() ); assert( w_diagonal_stiff_vector.isnonnegative() );
                assert( w_diagonal_stiff_pseudo.isfinite() ); assert( w_diagonal_stiff_pseudo.isnonnegative() );
                
                max_s_scalar_mass[r-r_min] = s_diagonal_mass_scalar.maxnorm();
                max_s_vector_mass[r-r_min] = s_diagonal_mass_vector.maxnorm();
                max_s_pseudo_mass[r-r_min] = s_diagonal_mass_pseudo.maxnorm();
                max_s_volume_mass[r-r_min] = s_diagonal_mass_volume.maxnorm();
                
                min_s_scalar_mass[r-r_min] = s_diagonal_mass_scalar.min();
                min_s_vector_mass[r-r_min] = s_diagonal_mass_vector.min();
                min_s_pseudo_mass[r-r_min] = s_diagonal_mass_pseudo.min();
                min_s_volume_mass[r-r_min] = s_diagonal_mass_volume.min();

                max_w_scalar_mass[r-r_min] = w_diagonal_mass_scalar.maxnorm();
                max_w_vector_mass[r-r_min] = w_diagonal_mass_vector.maxnorm();
                max_w_pseudo_mass[r-r_min] = w_diagonal_mass_pseudo.maxnorm();
                max_w_volume_mass[r-r_min] = w_diagonal_mass_volume.maxnorm();
                
                min_w_scalar_mass[r-r_min] = w_diagonal_mass_scalar.min();
                min_w_vector_mass[r-r_min] = w_diagonal_mass_vector.min();
                min_w_pseudo_mass[r-r_min] = w_diagonal_mass_pseudo.min();
                min_w_volume_mass[r-r_min] = w_diagonal_mass_volume.min();

                const Float threshold = shift + desired_closeness;

                LOG << max_w_scalar_mass[r-r_min] << space << PowerMethod( w_mass_scalar, repetitions ) << nl;

                max_s_scalar_stiff[r-r_min] = s_diagonal_stiff_scalar.maxnorm() - shift;
                max_s_vector_stiff[r-r_min] = s_diagonal_stiff_vector.maxnorm() - shift;
                max_s_pseudo_stiff[r-r_min] = s_diagonal_stiff_pseudo.maxnorm() - shift;
                
                min_s_scalar_stiff[r-r_min] = lowest( s_diagonal_stiff_scalar, threshold ) - shift;
                min_s_vector_stiff[r-r_min] = lowest( s_diagonal_stiff_vector, threshold ) - shift;
                min_s_pseudo_stiff[r-r_min] = lowest( s_diagonal_stiff_pseudo, threshold ) - shift;
                
                max_w_scalar_stiff[r-r_min] = w_diagonal_stiff_scalar.maxnorm() - shift;
                max_w_vector_stiff[r-r_min] = w_diagonal_stiff_vector.maxnorm() - shift;
                max_w_pseudo_stiff[r-r_min] = w_diagonal_stiff_pseudo.maxnorm() - shift;
                
                min_w_scalar_stiff[r-r_min] = lowest( w_diagonal_stiff_scalar, threshold ) - shift;
                min_w_vector_stiff[r-r_min] = lowest( w_diagonal_stiff_vector, threshold ) - shift;
                min_w_pseudo_stiff[r-r_min] = lowest( w_diagonal_stiff_pseudo, threshold ) - shift;
                
            }
            
        } 
    
        LOG << "\nSullivan mass matrix (max/min/cond)\n" << nl;
    
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e %.15e\n", r,
                max_s_scalar_mass[r-r_min], max_s_vector_mass[r-r_min], max_s_pseudo_mass[r-r_min], max_s_volume_mass[r-r_min] );
        
        LOG << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e %.15e\n", r,
                min_s_scalar_mass[r-r_min], min_s_vector_mass[r-r_min], min_s_pseudo_mass[r-r_min], min_s_volume_mass[r-r_min] );
        
        LOG << nl << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e %.15e\n", r,
                max_s_scalar_mass[r-r_min]/min_s_scalar_mass[r-r_min],
                max_s_vector_mass[r-r_min]/min_s_vector_mass[r-r_min],
                max_s_pseudo_mass[r-r_min]/min_s_pseudo_mass[r-r_min],
                max_s_volume_mass[r-r_min]/min_s_volume_mass[r-r_min]
                );
        
        LOG << "\nWhitney mass matrix (max/min/cond)\n" << nl;
    
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e %.15e\n", r,
                max_w_scalar_mass[r-r_min], max_w_vector_mass[r-r_min], max_w_pseudo_mass[r-r_min], max_w_volume_mass[r-r_min] );
        
        LOG << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e %.15e\n", r,
                min_w_scalar_mass[r-r_min], min_w_vector_mass[r-r_min], min_w_pseudo_mass[r-r_min], min_w_volume_mass[r-r_min] );
        
        LOG << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e %.15e\n", r,
                max_w_scalar_mass[r-r_min]/min_w_scalar_mass[r-r_min],
                max_w_vector_mass[r-r_min]/min_w_vector_mass[r-r_min],
                max_w_pseudo_mass[r-r_min]/min_w_pseudo_mass[r-r_min],
                max_w_volume_mass[r-r_min]/min_w_volume_mass[r-r_min]
                );
        
        LOG << nl;
        
        LOG << "\nSullivan stiff matrix (max/min/cond)\n" << nl;
    
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e\n", r,
                max_s_scalar_stiff[r-r_min], max_s_vector_stiff[r-r_min], max_s_pseudo_stiff[r-r_min] );
        
        LOG << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e\n", r,
                min_s_scalar_stiff[r-r_min], min_s_vector_stiff[r-r_min], min_s_pseudo_stiff[r-r_min] );
        
        LOG << nl << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e\n", r,
                max_s_scalar_stiff[r-r_min]/min_s_scalar_stiff[r-r_min],
                max_s_vector_stiff[r-r_min]/min_s_vector_stiff[r-r_min],
                max_s_pseudo_stiff[r-r_min]/min_s_pseudo_stiff[r-r_min]
                );
        
        LOG << "\nWhitney stiff matrix (max/min/cond)\n" << nl;
    
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e\n", r,
                max_w_scalar_stiff[r-r_min], max_w_vector_stiff[r-r_min], max_w_pseudo_stiff[r-r_min] );
        
        LOG << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e\n", r,
                min_w_scalar_stiff[r-r_min], min_w_vector_stiff[r-r_min], min_w_pseudo_stiff[r-r_min] );
        
        LOG << nl;
        
        for( int r = r_min; r <= r_max; r++ ) 
            LOGPRINTF("%i:\t%.15e %.15e %.15e\n", r,
                max_w_scalar_stiff[r-r_min]/min_w_scalar_stiff[r-r_min],
                max_w_vector_stiff[r-r_min]/min_w_vector_stiff[r-r_min],
                max_w_pseudo_stiff[r-r_min]/min_w_pseudo_stiff[r-r_min]
                );
        
        LOG << nl;
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
