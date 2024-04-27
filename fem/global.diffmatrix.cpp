#include <vector>

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "global.diffmatrix.hpp"

SparseMatrix FEECBrokenDiffMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 1 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <  n );
    
    
    // generate an empty local sample matrix
    
    const std::vector<MultiIndex> multis_dest = generateMultiIndices( IndexRange( 0, n ), r-1 );
    const std::vector<MultiIndex> multis_src  = generateMultiIndices( IndexRange( 0, n ), r );
    
    const std::vector<IndexMap>   sigmas_dest = generateSigmas( IndexRange( 1, k+1 ), IndexRange( 0, n ) );
    const std::vector<IndexMap>   sigmas_src  = generateSigmas( IndexRange( 1, k   ), IndexRange( 0, n ) );
    
    const int localdim_out = multis_dest.size() * sigmas_dest.size();
    const int localdim_in  = multis_src.size()  * sigmas_src.size();
    
    std::vector<SparseMatrix::MatrixEntry> localmatrixentries;
    
    assert( multis_dest.size() == binomial_integer( n+r-1, n   ) );
    assert( multis_src.size()  == binomial_integer( n+r,   n   ) );
    assert( sigmas_dest.size() == binomial_integer( n+1,   k+1 ) );
    assert( sigmas_src.size()  == binomial_integer( n+1,   k   ) );
    
    for( int src_poly_index = 0; src_poly_index < multis_src.size(); src_poly_index++ )
    for( int src_form_index = 0; src_form_index < sigmas_src.size(); src_form_index++ )
    for( int p = 0; p <= n; p++ )
    {
        
        const MultiIndex& src_poly = multis_src[src_poly_index];
        const IndexMap&   src_form = sigmas_src[src_form_index];
        
        if( src_poly[p] == 0 or src_form.has_value_in_range(p) )
            continue;
        
        MultiIndex new_poly = src_poly - p;
        IndexMap   new_form = expand_one( src_form, p );
        
        int new_poly_index = find_index( multis_dest, new_poly );
        int new_form_index = find_index( sigmas_dest, new_form );
        
        assert( new_form.getSourceRange().min() == 1 );
        int signum = sign_power( new_form.preimageof(p) - 1 );
        
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = new_poly_index * sigmas_dest.size() + new_form_index;
        entry.column = src_poly_index * sigmas_src.size()  + src_form_index;
        entry.value  = src_poly[p] * signum;
        
        assert( entry.row >= 0 && entry.row < localdim_out );
        assert( entry.column >= 0 && entry.column < localdim_in );
        
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



