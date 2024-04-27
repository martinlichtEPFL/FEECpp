
#include <vector>

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/simpleoperators.hpp"
#include "../mesh/mesh.hpp"

#include "indexfunctions.hpp"

#include "../fem/global.flags.hpp"



DiagonalOperator FEECSullivanFlagMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r >= 1 or k == n ); // NB: we allow r == 0 in case of volume forms 

    // generate the list of sigmas and multiindices for each dimension 
    
    std::vector< std::vector< std::pair<MultiIndex,IndexMap> > > lists_of_sullivan_indices( n+1 );
    for( auto d : IndexRange(0,n) )
        lists_of_sullivan_indices[d] = ListOfSullivanIndices( d, k, r );

    // count mesh properties 
    
    std::vector< int > num_faces( n+1 );
    for( auto d : IndexRange(0,n) )
        num_faces[d] = mesh.count_simplices( d );

    // for auxiliary purposes, compute the column offsets 
    
    std::vector<int> column_offset( n+1, 0 );
    for( int d = 1; d <= n; d++ )
        column_offset[d] = column_offset[d-1] + mesh.count_simplices(d-1) * lists_of_sullivan_indices[d-1].size();
    
    // dimensions of the matrix 
    
    int dim_continuous  = 0;
    for( auto d : IndexRange(0,n) )
        dim_continuous += num_faces[d] * lists_of_sullivan_indices[d].size();
    
    
    // create that matrix 
    
    DiagonalOperator ret( dim_continuous, 0. );
    
    for( int d = 0; d <= n; d++ ) // go over all the subsimplex dimensions
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif 
    for( int f = 0; f < num_faces[d]; f++ ) // go over all the d dimensional subsimplices 
    for( int index_alphasigma = 0; index_alphasigma < lists_of_sullivan_indices[d].size(); index_alphasigma++ )   // go over the corresponding alpha/sigma pairs
    {
        
        int continuous_index = column_offset[d] + f * lists_of_sullivan_indices[d].size() + index_alphasigma;

        assert( 0 <= continuous_index and continuous_index < dim_continuous  );

        ret.getdiagonal().setentry( continuous_index, mesh.get_flag( d, f ) == SimplexFlag::SimplexFlagDirichlet ? 0. : 1. );
        
    }
    
    return ret;
}





DiagonalOperator FEECWhitneyFlagMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r >= 1 );
    
    // generate the list of rhos and multiindices for each dimension
    
    std::vector< std::vector< std::pair<MultiIndex,IndexMap> > > lists_of_Whitney_indices( n+1 );
    for( auto d : IndexRange(0,n) )
        lists_of_Whitney_indices[d] = ListOfWhitneyIndices( d, k, r );

    // count mesh properties 
    
    std::vector< int > num_faces( n+1 );
    for( auto d : IndexRange(0,n) )
        num_faces[d] = mesh.count_simplices( d );

    // dimensions of the matrix 
    
    int dim_continuous  = 0;
    for( auto d : IndexRange(0,n) )
        dim_continuous += num_faces[d] * lists_of_Whitney_indices[d].size();
    
    // for auxiliary purposes, compute the column offsets 
    
    std::vector<int> column_offset( n+1, 0 );
    for( int d = 1; d <= n; d++ )
        column_offset[d] = column_offset[d-1] + mesh.count_simplices(d-1) * lists_of_Whitney_indices[d-1].size();

    
    // create the matrix 

    DiagonalOperator ret( dim_continuous, 0. );
    
    for( int d = 0; d <= n; d++  ) // go over all the subsimplex dimensions
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif 
    for( int f = 0; f < num_faces[d]; f++  ) // go over all the d dimensional simplices 
    for( int index_alpharho = 0; index_alpharho < lists_of_Whitney_indices[d].size(); index_alpharho++ )
    {
        
        int continuous_index = column_offset[d] + f * lists_of_Whitney_indices[d].size() + index_alpharho;

        assert( 0 <= continuous_index and continuous_index < dim_continuous  );

        ret.getdiagonal().setentry( continuous_index, mesh.get_flag( d, f ) == SimplexFlag::SimplexFlagDirichlet ? 0. : 1. );

    }
    
    return ret;
}

