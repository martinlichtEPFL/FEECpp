
#include <vector>

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/linearoperator.hpp"
#include "../mesh/mesh.hpp"

#include "indexfunctions.hpp"

#include "../fem/global.afwincl.hpp"


SparseMatrix FEECAFWInclusionMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    // check whether the parameters are right 
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    assert( r >= 1 or k == n ); // NB: we allow r == 0 in case of volume forms 
    
    // generate the list of rhos and multiindices for each dimension
    
    std::vector< std::vector< std::pair<MultiIndex,IndexMap> > > lists_of_Whitney_indices( n+1 );
    for( auto d : IndexRange(0,n) )
        lists_of_Whitney_indices[d] = ListOfWhitneyIndices( d, k, r );

    const auto multis_dst = generateMultiIndices( IndexRange(0,n), r );
    const auto sigmas_dst = generateSigmas( IndexRange(1,k), IndexRange(0,n) );
    
    
    // count mesh properties 
    
    const int num_volumes = mesh.count_simplices( n );
    
    std::vector< int > num_faces( n+1 );
    for( auto d : IndexRange(0,n) )
        num_faces[d] = mesh.count_simplices( d );

    
    // dimensions of the matrix 
    
    const int dim_broken = num_volumes * binomial_integer( n+1, k ) * binomial_integer( n+r, r );
    
    int dim_continuous  = 0;
    for( auto d : IndexRange(0,n) )
        dim_continuous += num_faces[d] * lists_of_Whitney_indices[d].size();
    
    int num_entries = 0;
    for( auto d : IndexRange(0,n) )
        num_entries += num_volumes * binomial_integer(n+1,d+1) * lists_of_Whitney_indices[d].size() * (k+1);
    
    // create that matrix 
    
    SparseMatrix ret( dim_broken, dim_continuous, num_entries );
    
    
    // for auxiliary purposes, create the index inclusion maps 
    // from subsimplices into the supersimplices
    
    std::vector< std::vector<IndexMap> > subsimplex_inclusions( n+1 );
    for( auto d : IndexRange(0,n) )
        subsimplex_inclusions[d] = generateSigmas( IndexRange(0,d), IndexRange(0,n) );
    

    // for auxiliary purposes, compute the column offsets 
    std::vector<int> column_offset( n+1, 0 );
    for( int d = 1; d <= n; d++ )
        column_offset[d] = column_offset[d-1] + mesh.count_simplices(d-1) * lists_of_Whitney_indices[d-1].size();

    // for auxiliary purposes, compute the entries offsets 
    std::vector<int> entries_offset( n+1, 0 );
    for( int d = 1; d <= n; d++ )
        entries_offset[d] = entries_offset[d-1] + binomial_integer(n+1,d) * lists_of_Whitney_indices[d-1].size() * (k+1) * num_volumes;


    for( int d  = 0; d  <= n;                          d++  )    // go over all the subsimplex dimensions
    for( int fi = 0; fi  < binomial_integer(n+1,d+1); fi++  )    // go over all the d dimensional subsimplices 
    // for( const auto& alpharho : lists_of_Whitney_indices[d] )    // go over the corresponding alpha/rho pairs
    for( int index_alpharho = 0; index_alpharho < lists_of_Whitney_indices[d].size(); index_alpharho++ )
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif 
    for( int s  = 0; s   < num_volumes;                s++  )    // go over all the volumes 
    {
        
        // Find indices of things and prepare auxiliary variables 
        
        const auto& alpharho = lists_of_Whitney_indices[d][index_alpharho];

        const MultiIndex& alpha = alpharho.first;
        
        const IndexMap&     rho = alpharho.second;
        
        // const int index_alpharho = 
        //     std::find( lists_of_Whitney_indices[d].begin(), lists_of_Whitney_indices[d].end(), alpharho )
        //     - 
        //     lists_of_Whitney_indices[d].begin();
        assert( 0 <= index_alpharho && index_alpharho < lists_of_Whitney_indices[d].size() );
        assert( lists_of_Whitney_indices[d][index_alpharho] == alpharho );
        
        const int index_fi = mesh.get_subsimplex( n, d, s, fi );
        assert( 0 <= index_fi && index_fi < mesh.count_simplices(d) );
        
        
        // get inclusion index map
                
        IndexMap volume_vertices = mesh.getsubsimplices( n, 0, s        );
        IndexMap face_vertices   = mesh.getsubsimplices( d, 0, index_fi );
        int inclusion_index = 0;
        for( ; inclusion_index < subsimplex_inclusions[d].size(); inclusion_index++ )
            if( face_vertices == volume_vertices * subsimplex_inclusions[d][inclusion_index] )
                break;
        assert( inclusion_index < subsimplex_inclusions[d].size() );
        assert( face_vertices == volume_vertices * subsimplex_inclusions[d][inclusion_index] );
        const IndexMap inclusion = subsimplex_inclusions[d][inclusion_index];
        

        // create the matrix entries 

        for( int j = 0; j <= k; j++ )
        {

            // const MultiIndex alpha_effective = alpha;
            // alpha_effective.add( rho[j] );
            const MultiIndex alpha_effective = [&,alpha,rho,j]()->MultiIndex{
                auto temp = alpha;
                temp.add( rho[j] );
                return temp;
            }();

            // IndexMap sigma_effective = rho.skip(j).shiftup();
            const IndexMap sigma_effective = [&,rho,j]()->IndexMap{
                auto values = rho.getvalues();
                values.erase( values.begin() + j );
                return IndexMap( IndexRange(1,k), IndexRange(0,d), values );
            }();
            
            
            // create actual multiindices 
            
            const MultiIndex alpha_vol = MultiIndex( IndexRange(0,n), [&alpha_effective,&inclusion]( int p ) -> int {
                                                assert( inclusion.getTargetRange().contains(p) ); 
                                                if( inclusion.has_value_in_range(p) )
                                                    return alpha_effective.at( inclusion.preimageof(p) );
                                                else
                                                    return 0;
                                            } );
            
            const IndexMap   sigma_vol = inclusion * sigma_effective;
            
            
            // find the indices of those extended functions 
            
            const auto index_alpha_vol = std::find( multis_dst.begin(), multis_dst.end(), alpha_vol ) - multis_dst.begin();
            
            const auto index_sigma_vol = std::find( sigmas_dst.begin(), sigmas_dst.end(), sigma_vol ) - sigmas_dst.begin();

//             if( 0 <= index_sigma_vol and index_sigma_vol < sigmas_dst.size() );
//             else std::cout << rho << tab << j << tab << sigma_effective << sigma_vol << nl;
//             if( 0 <= index_alpha_vol and index_alpha_vol < multis_dst.size() );
//             else std::cout << alpha_vol << nl;
            


            assert( 0 <= index_alpha_vol and index_alpha_vol < multis_dst.size() );
            assert( 0 <= index_sigma_vol and index_sigma_vol < sigmas_dst.size() );
            

            // enter the values of the data structure 

            int broken_index = s * binomial_integer(n+r,r) * binomial_integer(n+1,k) 
                           +
                           index_alpha_vol * binomial_integer(n+1,k)
                           + 
                           index_sigma_vol;

            int continuous_index = column_offset[d] // sum_int( d-1, [&mesh,&lists_of_Whitney_indices](int i) -> int { return mesh.count_simplices(i) * lists_of_Whitney_indices[i].size(); } )
                           + 
                           index_fi * lists_of_Whitney_indices[d].size()
                           +
                           index_alpharho
                           ;

            Float value  = sign_power(j);

            int index_of_entry = entries_offset[d] // sum_int( d-1, [ &lists_of_Whitney_indices, n ](int c) -> int { return binomial_integer(n+1,c+1) * lists_of_Whitney_indices[c].size(); } ) * (k+1) * num_volumes
                                 +
                                 fi * lists_of_Whitney_indices[d].size() * (k+1) * num_volumes
                                 +
                                 index_alpharho * (k+1) * num_volumes
                                 +
                                 j * num_volumes
                                 +
                                 s; 
                                
            assert( broken_index     < dim_broken );
            assert( continuous_index < dim_continuous  );
            assert( index_of_entry   < num_entries );
            assert( broken_index >= 0 && continuous_index >= 0 && index_of_entry >= 0 );
                
            // set up the actual entry
            
            SparseMatrix::MatrixEntry entry;
            
            entry.row    = broken_index;
            entry.column = continuous_index;
            entry.value  = value;
            
            if( mesh.get_flag( d, index_fi ) == SimplexFlag::SimplexFlagDirichlet )
                entry.value = 0.;
            
            ret.setentry( index_of_entry, entry );

        }
        
        
    }
    
    return ret;
}

