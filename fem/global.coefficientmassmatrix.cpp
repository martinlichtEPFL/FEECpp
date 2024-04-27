
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/factorization.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"
#include "../fem/utilities.hpp"

#include "../fem/polynomialmassmatrix.hpp"

#include "../fem/global.coefficientmassmatrix.hpp"



SparseMatrix FEECBrokenCoefficientMassMatrix( const Mesh& mesh, int n, int k, int r,
                                              int w, const std::function<DenseMatrix(const FloatVector&)>& generator 
) {
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( binomial_integer( n+r, n ) == binomial_integer( n+r, r ) );
    assert( w >= 0 );
    
    // Dimensions of the output matrix and number of entries 
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int dim_in      = num_simplices * localdim;
    const int dim_out     = num_simplices * localdim;
    const int num_entries = num_simplices * localdim * localdim;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    
    // assemble algebraic auxiliary material
    // - lagrange points in barycentric coordinates 
    // - coefficients of Lagrange polynomials
    // mass matrices 

    const auto lpbcs = InterpolationPointsInBarycentricCoordinates( n, w );

    const auto lpcoeff = Inverse( PointValuesOfMonomials( w, lpbcs ) );
    
    const auto polymassmatrix_per_point = polynomialmassmatrices_per_lagrangepoint( n, r, w );
    
    
    // loop over the simplices and compute the mass matrices

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        // assemble some data for the element 
        // - measure 
        // - barycentric coordinates 
        // - lagrange points 
        
        const Float measure     = mesh.getMeasure( n, s );
        assert( measure >= 0. );

        const DenseMatrix GM    = mesh.getGradientMatrix( n, s );
        const DenseMatrix extGM = SubdeterminantMatrix( GM, k );

        const auto vertex_coordinates = mesh.getVertexCoordinateMatrix( n, s );
        const auto lpeucl             = vertex_coordinates * lpbcs;

        // compute the mass matrix contribution 
        // for each lagrange point 
        
        DenseMatrix full_element_matrix( localdim, localdim, 0. );

        for( int p = 0; p < polymassmatrix_per_point.size(); p++ )
        {
            const auto& polyMM = polymassmatrix_per_point[p];
            
            const DenseMatrix matrix_at_point = generator( lpeucl.getcolumn(p) );

            const auto formMM = Transpose(extGM) * matrix_at_point * extGM;
            // const auto formMM = MatrixTripleMult( matrix_at_point, extGM );

            // DenseMatrix GPM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, s ), k );
            // assert( ( GPM - formMM ).is_numerically_small() );

            if( w == 0 ) assert( ( polyMM - polynomialmassmatrix(n,r) ).is_numerically_small() );

            //auto fullMM = measure * MatrixTensorProduct( polyMM, formMM );
            auto fullMM = MatrixTensorProduct( polyMM, formMM );
            fullMM *= measure;

            full_element_matrix += fullMM;
        }
        
        // DONE ... now list everything.

        for( int row = 0; row < localdim; row++ )
        for( int col = 0; col < localdim; col++ )
        {
            int index_of_entry = s * localdim * localdim + row * localdim + col;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * localdim + row;
            entry.column = s * localdim + col;
            entry.value  = full_element_matrix( row, col );
            
            ret.setentry( index_of_entry, entry );
        }

    }

    LOG << "Finished Sparse Matrix entries: " << num_entries << "\n";
    
    return ret;
}




