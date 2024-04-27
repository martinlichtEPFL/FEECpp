
#include <vector>

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"
#include "../fem/utilities.hpp"

#include "../fem/global.interpol.hpp"


SparseMatrix FEECBrokenInterpolationMatrix( const Mesh& mesh, int n, int k, int r, int r_plus )
{
    
    // check whether the parameters are right 
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r_plus >= 0 );
    
    const std::vector<MultiIndex> multis_low  = generateMultiIndices( IndexRange( 0, n ), r          );
    const std::vector<MultiIndex> multis_high = generateMultiIndices( IndexRange( 0, n ), r + r_plus );
    
    assert( multis_low.size()  == binomial_integer( n + r         , n ) );
    assert( multis_high.size() == binomial_integer( n + r + r_plus, n ) );
    
    const std::vector<IndexMap> sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n ) );
    
    assert( sigmas.size() == binomial_integer( n+1, k ) );


    // check whether the parameters are right 

    const DenseMatrix bcs = InterpolationPointsInBarycentricCoordinates( n, r );

    const DenseMatrix pvom_plus = PointValuesOfMonomials( r + r_plus, bcs );

    const DenseMatrix pvom      = PointValuesOfMonomials( r         , bcs );

    assert( pvom.issquare() );

    const DenseMatrix lpc = LagrangePolynomialCoefficients( n, r );
    const DenseMatrix lpc_inv = Inverse(lpc);

    const DenseMatrix localmatrix = MatrixTensorProduct( Inverse( pvom ) * pvom_plus, IdentityMatrix( sigmas.size() ) );


    
    

    
    const int localdim_in  = multis_high.size() * sigmas.size();
    const int localdim_out = multis_low.size()  * sigmas.size();
    
    std::vector<SparseMatrix::MatrixEntry> localmatrixentries;
    
    for( int  low_poly_index = 0; low_poly_index  <  multis_low.size();  low_poly_index++ )
    for( int high_poly_index = 0; high_poly_index < multis_high.size(); high_poly_index++ )
    for( int      form_index = 0;      form_index <      sigmas.size();      form_index++ )
    {
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = low_poly_index  * sigmas.size() + form_index;
        entry.column = high_poly_index * sigmas.size() + form_index;
        entry.value  = localmatrix( entry.row, entry.column );
        
        assert( entry.row    >= 0 && entry.row    < localdim_out );
        assert( entry.column >= 0 && entry.column < localdim_in  );
        
        localmatrixentries.push_back( entry );
        
    }
    
    // Auxiliary calculations and preparations
    
    // Finished generating local matrix
    
    const int num_simplices = mesh.count_simplices( n );
        
    int noe = localmatrixentries.size();
    
    const int dim_out = num_simplices * localdim_out;
    const int dim_in  = num_simplices * localdim_in;
    
    const int num_entries = num_simplices * noe;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices;   s++ )
    for( int i = 0; i < noe; i++ )
    {
        
        int index_of_entry = s * noe + i; 
            
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = s * localdim_out + localmatrixentries[i].row;
        entry.column = s * localdim_in  + localmatrixentries[i].column;
        entry.value  = localmatrixentries[i].value;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}

