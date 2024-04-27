#include <algorithm>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/utilities.hpp"


// 
// Generates a matrix whose columns are the barycentric coordinates
// of the interpolation points over an n dimensional simplex 
// such that polynomials of degree r can be reconstructed.
// 
// The barycentric coordinates are chosen such that 
// the points are located in the interior of the simplex.
// (See parameter delta)
// 
// Size of returned matrix:
// [n+1] x [ n+r choose r ]

DenseMatrix InterpolationPointsInBarycentricCoordinates( int n, int r )
{
    assert( 0 <= n && 0 <= r );
    
    const auto multi_indices = generateMultiIndices( IndexRange(0,n), r );
    
    DenseMatrix ret( n+1, SIZECAST( multi_indices.size() ) );
    
    assert( multi_indices.size() == binomial_integer( n+r, r ) );
    for( const auto& mi : multi_indices )
        assert( mi.absolute() == r );
    assert( ret.getdimout() == n+1 );
    assert( ret.getdimin() == multi_indices.size() );
    
    const Float delta = 0.1;
    
    assert( delta > 0 );
    
    for( int i = 0; i < multi_indices.size(); i++ )
        ret.setcolumn( i, FloatVector( multi_indices[i].getvalues() ).shift( delta ).scaleinverse( r + (n+1) * delta ) );
    // The scaling factor ensures that the entries define a convex combination 

    assert( ret.isfinite() );
    for( int i = 0; i < ret.getdimin(); i++ ) {
        assert( ret.getcolumn(i).isnonnegative() );
        assert( ret.getcolumn(i).sumnorm() > 0.9999 && ret.getcolumn(i).sumnorm() < 1.0001 );
    }
    
    return ret;
}




// 
// Suppose that the columns of bcs are evaluation points
// (in barycentric coordinates) over a d simplex.
// 
// The output matrix is as follows:
// - row indices correspond to the evaluation points 
// - column indices correspond to the multiindices (standard Lagrange basis of degree r)
// - the entries are value of the corresponding polynomial at the corresponding point
// 
// Size of returned matrix:
// [ number of evaluation points ] x [ n+r choose r ]
// 
// NOTE: if M is the output matrix and x some coefficient vector for a polynomial,
//       then y = M * x is the values of that polynomial at the Lagrange points
//       If M is invertible, then x = inv(M) y provides the monomial coefficients 
//       for a given vector y of point values.

DenseMatrix PointValuesOfMonomials( int r, const DenseMatrix& bcs )
{
    assert( 0 <= r );
    assert( 0 < bcs.getdimout() );
    
    const int dim = bcs.getdimout()-1;
    
    const auto mis = generateMultiIndices( IndexRange(0,dim), r );
    
    const int number_of_evaluation_points = bcs.getdimin();
    
    const int number_of_polynomials = mis.size();
    
    assert( mis.size() == binomial_integer( dim+r, r ) );
    
    DenseMatrix ret( number_of_evaluation_points, number_of_polynomials );
    
    for( int row = 0; row < number_of_evaluation_points; row++ ) // row -> interpolation point 
    for( int col = 0; col < number_of_polynomials;       col++ ) // col -> barycentric poly 
    {
        ret(row,col) = 1.;
        for( int d = 0; d <= dim; d++ )
            if( mis[col][d] != 0 )
                ret(row,col) *= power_numerical( bcs(d,row), mis[col][d] );
    }
    
    assert( ret.isfinite() );
    
    return ret;
    
}




// 
// Suppose that the columns of bcs are evaluation points
// (in barycentric coordinates) over a d simplex.
// 
// The output matrix is as follows:
// - row indices correspond to the multiindices (standard Lagrange basis of degree r)
// - column indices correspond to the evaluation points 
// - each entry is the coefficient of the corresponding multindex 
//   for the Lagrange polynomial associated to the evaluation point 
// 
// Size of returned matrix:
// [ n+r choose r ] x [ n+r choose r ]

DenseMatrix LagrangePolynomialCoefficients( int n, int r )
{
    auto lagrangepoints_baryc = InterpolationPointsInBarycentricCoordinates( n, r );

    return Inverse( PointValuesOfMonomials( r, lagrangepoints_baryc ) );
}




// 
// The input matrix J is a Jacobian of the transformation which transforms 
// the [dimin] reference simplex onto a simplex 
// within [dimout] dimensional ambient space 
// 
// We transpose the matrix and prepend a zero row.
// 
// If we multiply a vector in standard basis coordinates 
// we get the same vector in barycentric coordinates, 
// using only the gradients 1 through n. 
// Gradient 0 is redundant and thus the zero-th row is 0.
// 
// Size of output matrix:
// [ J.dimin()+1 x J.dimout() ]
// 
// TODO: Remove this function 
// 

// DenseMatrix BarycentricProjectionMatrix( const DenseMatrix& J )
// {
//     assert( J.getdimout() >= J.getdimin() );
    
//     const auto F = Transpose(J);
    
//     DenseMatrix ret( J.getdimin()+1, J.getdimout(), 0.0 );
    
//     for( int r = 0; r < J.getdimin();  r++ )
//     for( int c = 0; c < J.getdimout(); c++ )
//         ret( r+1, c ) = F(r,c);
    
//     assert( ret.isfinite() );
    
//     return ret;
// }




// 
// evaluate a field at given physical points 
// and collect its values 
// 
// Given a k-differential form in dim dimensional ambient space 
// and lagrangepoints_eucl being a matrix of size [outerdimension] x [var1]
// we let the output be as follows.
// 
// Here, [var1] is the number of evaluation points.
// 
// The output is a matrix of size [dim choose k] x [var1]
// whose column are the evaluations of the field 
// at each of the evaluation points, given in Euclidean coordinates.
// 

DenseMatrix EvaluateField( 
            int outerdim, int k,  
            const DenseMatrix& lagrangepoints_eucl, 
            const std::function< FloatVector( const FloatVector& ) >& field
            )
{
    
    assert( 0 <= outerdim );
    assert( 0 <= k && k <= outerdim );
    assert( lagrangepoints_eucl.getdimout() == outerdim );
    assert( lagrangepoints_eucl.isfinite() );
    
    const auto fielddim = binomial_integer(outerdim,k);
    
    const auto number_of_evaluation_points = lagrangepoints_eucl.getdimin();
    
    DenseMatrix ret( fielddim, number_of_evaluation_points );
    
    for( int p = 0; p < number_of_evaluation_points; p++ )
    {
        const auto evaluation_point = lagrangepoints_eucl.getcolumn(p);

        assert( evaluation_point.isfinite() );
        
        const auto value = field( evaluation_point );
        
        Assert( value.isfinite(), evaluation_point, value );
        
        assert( value.getdimension() == binomial_integer( outerdim, k ) );
        
        ret.setcolumn( p, value );
    }
    
    return ret;
    
}





// TODO:
// 
// Given a function in ambient space, giving forms in Euclidean coordinates
// 
// 1. Evaluate the form at the Lagrange points 
// 
// 2. Get the coefficients of the barycentric polynomials
// 
// 3. Transform the Euclidean components to their barycentric projections
// 
// The output is a vector whose coefficients represent the volume-wise interpolations 
// into the local Sullivan space of k-forms with polynomial degree r 
// 

FloatVector Interpolation( 
            const Mesh& m, 
            int dim, int k, int r, 
            const std::function< FloatVector( const FloatVector& ) >& field
            )
{
    
    assert( 0 <= dim && dim <= m.getinnerdimension() );
    assert( 0 <= k && k <= dim );
    assert( 0 <= r );
    
    FloatVector ret( m.count_simplices(dim) * SullivanSpanSize(dim,k,r) );
    
    const int outerdim = m.getouterdimension();
    
    // lagrangepoints_baryc are the barycentric coordinates to interpolate 
    // degree r polynomials over a dim simplex 
    // [dim+1] x [number of lagrange points] 
    // 
    // M is the function values of the Lagrange polynomials at those points 
    // [n+r choose r] x [number of lagrange points] (actually, square)
    // 
    // Minv is the inverse of M.
    // [same size as M]
    
    const auto lagrangepoints_baryc = InterpolationPointsInBarycentricCoordinates( dim, r );
    
    const auto M = PointValuesOfMonomials( r, lagrangepoints_baryc );
        
    const auto Minv = Inverse( M );
    
    assert( lagrangepoints_baryc.getdimin()  == SullivanSpanSize(dim,0,r) );
    assert( lagrangepoints_baryc.getdimout() == dim+1 );
    assert( M.getdimout()    == SullivanSpanSize(dim,0,r) );
    assert( M.getdimin()     == SullivanSpanSize(dim,0,r) );
    assert( M.isfinite()    );
    assert( Minv.isfinite() );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < m.count_simplices(dim); s++ )
    {
        
        // We obtain the coordinate matrix (outerdim) x (dim+1) 
        // of the dim-dimensional simplex with index s.
        // Each column contains the physical coordinates of the vertices
        // 
        // lagrangepoints_eucl are the physical coordinates of the interpolation points 
        // (outerdim) x (number of Lagrange points)
        

        const auto euclidean_coords = m.getVertexCoordinateMatrix( dim, s );
        
        const auto lagrangepoints_eucl = euclidean_coords * lagrangepoints_baryc;
        
        
        // Jac is a matrix of size [outerdimension] x [dim]
        // which is the jacobian of the transformation mapping of simplex s
        // 
        // bpm is a matrix of size [dim+1]x[outerdimension]
        // that is the barycentric projection matrix (see paper)
        // 
        // P is the matrix of subdeterminants of size k of that matrix 
        
        
        const auto bpm = m.getBarycentricProjectionMatrix( dim, s );

        const auto P = SubdeterminantMatrix( bpm, k );
        
        
        // TODO: complete remove the following two lines 
        // const auto Jac = m.getTransformationJacobian( dim, s );
        // const auto bpm_ = BarycentricProjectionMatrix( Jac );
        // assert( ( bpm - bpm_ ).norm() == 0. );
        
        // The InterpolationMatrix is the tensor product matrix of Minv and P.
        // 
        // Evaluations is a matrix of size [dim choose k] x [lagrangepoints_eucl.dimin]
        // whose columsn are the outputs of the field at the Lagrange points
        // 
        // localResults contains the coeffecients of the interpolation
        // in the canonical basis (Lagrange polynomials x barycentric gradients)
        
        const auto InterpolationMatrix = MatrixTensorProduct( Minv, P );
        
        const auto Evaluations      = EvaluateField( outerdim, k, lagrangepoints_eucl, field );
        const auto EvaluationVector = Evaluations.flattencolumns();
        
        const auto localResult = InterpolationMatrix * EvaluationVector;

        #ifndef NDEBUG

        assert( lagrangepoints_eucl.isfinite() );
        assert( P.isfinite() );
        assert( InterpolationMatrix.isfinite() );
        assert( Evaluations.isfinite()      );
        assert( EvaluationVector.isfinite() );
        
        assert( localResult.isfinite() );

        if( k == 0 ) {
            
            assert( ( P - DenseMatrix( 1, 1, 1. ) ).iszero() );
            assert( ( InterpolationMatrix - Minv ).iszero() );
            assert( ( InterpolationMatrix - Minv ).is_numerically_small() );
            
            if( not ( M * localResult - EvaluationVector ).is_numerically_small() ) {
                LOG << M << nl << nl;
                LOG << Minv << nl << nl;
                LOG << InterpolationMatrix << nl << nl;
                LOG << M * Minv << nl << nl;
                LOG << M * InterpolationMatrix << nl << nl;
                LOG << M * localResult << nl << nl << EvaluationVector << nl << nl;
                LOG << M * localResult - EvaluationVector << nl << nl;
                LOG << ( M * localResult - EvaluationVector ).norm() << nl << nl;
            }
            assert( ( M * localResult - EvaluationVector ).is_numerically_small() );
        }
        #endif
        
        assert( localResult.getdimension() == SullivanSpanSize(dim,k,r) );
        
        ret.setslice( s * SullivanSpanSize(dim,k,r), localResult );
        
    }
    
    return ret;
}

// std::vector<DenseMatrix> Interpolation( 
//             const Mesh& m, 
//             int dim, int r, 
//             const std::function< DenseMatrix( const FloatVector& ) >& matrixfield
//             )
// {
    
//     assert( 0 <= dim && dim <= m.getinnerdimension() );
//     assert( 0 <= r );
    
//     // std::vector<DenseMatrix> ret( m.count_simplices(dim) * binomial_integer(dim+r,r) );
    
//     const int outerdim = m.getouterdimension();
    
    
//     const auto lagrangepoints_baryc = InterpolationPointsInBarycentricCoordinates( dim, r );
    
//     const auto M = PointValuesOfMonomials( r, lagrangepoints_baryc );
        
//     const auto Minv = Inverse( M );
    
//     // TODO
    
//     // return ret;
// }




