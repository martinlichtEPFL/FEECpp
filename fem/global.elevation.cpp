
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
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/global.elevation.hpp"


SparseMatrix FEECBrokenElevationMatrix( const Mesh& mesh, int n, int k, int r, int r_plus )
{
    
    // check whether the parameters are right 
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r_plus >= 0 );
    
    const std::vector<MultiIndex> multis_add  = generateMultiIndices( IndexRange( 0, n ), r_plus     );
    const std::vector<MultiIndex> multis_low  = generateMultiIndices( IndexRange( 0, n ), r          );
    const std::vector<MultiIndex> multis_high = generateMultiIndices( IndexRange( 0, n ), r + r_plus );
    
    assert( multis_add.size()  == binomial_integer( n + r_plus    , n ) );
    assert( multis_low.size()  == binomial_integer( n + r,          n ) );
    assert( multis_high.size() == binomial_integer( n + r + r_plus, n ) );
    
    const std::vector<IndexMap> sigmas = generateSigmas( IndexRange( 1, k ), IndexRange( 0, n ) );
    
    assert( sigmas.size() == binomial_integer( n+1, k ) );
    
    const int localdim_in  = multis_low.size()  * sigmas.size();
    const int localdim_out = multis_high.size() * sigmas.size();
    
    std::vector<SparseMatrix::MatrixEntry> localmatrixentries;
    
    for( int low_poly_index = 0; low_poly_index < multis_low.size(); low_poly_index++ )
    for( int     form_index = 0;     form_index <     sigmas.size();     form_index++ )
    for( int add_poly_index = 0; add_poly_index < multis_add.size(); add_poly_index++ )
//     for( const MultiIndex& addendum : multis_add )
    {
        
        const MultiIndex& addendum = multis_add[add_poly_index];
        
        const MultiIndex& low_poly = multis_low[low_poly_index];
        
        MultiIndex high_poly = low_poly + addendum;
        
        int high_poly_index = find_index( multis_high, high_poly );
        
        assert( 0 <= high_poly_index and high_poly_index < multis_high.size() );
        
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = high_poly_index * sigmas.size() + form_index;
        entry.column = low_poly_index  * sigmas.size() + form_index;
        
        entry.value  = factorial_numerical( r_plus ) / addendum.factorial_numerical();
        
        assert( entry.row    >= 0 && entry.row    < localdim_out );
        assert( entry.column >= 0 && entry.column < localdim_in  );
        
        localmatrixentries.push_back( entry );
        
    }
    
// //     if( k==0 and r==0 and r_plus == 3 ){
// //         
// //         for( const auto& lme : localmatrixentries )
// //             LOG << lme.row << space << lme.column << space << lme.value;// << nl;
// //         
// //         exit(0);
// //     }
        
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

