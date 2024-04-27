
#include <algorithm>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/global.cechmatrix.hpp"

SparseMatrix FEECCechMassMatrix( const Mesh& mesh, int n, int k, int ss )
{
    
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    
    // Auxiliary calculations and preparations
    
    const int num_volumes = mesh.count_simplices( n );
    
    const int num_subcells = mesh.count_simplices( k );
    
    const int subcells_per_simplex = binomial_integer( n + 1, k + 1 );
    
    const int dim_in      = num_subcells;
    const int dim_out     = num_subcells;
    const int num_entries = num_volumes * subcells_per_simplex;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int volume = 0; volume < num_volumes;          volume++ )
    for( int cell   = 0; cell   < subcells_per_simplex; cell++   )
    {
        
        const Float measure = mesh.getMeasure( n, volume );
        assert( measure >= 0. );

        const int index = mesh.get_subsimplex( n, k, volume, cell );
        assert( 0 <= index && index < mesh.count_simplices(k) );
        
        Float diameter = 0.;

        if( k > 0 ) {
            diameter = mesh.getDiameter( k, index );
        } else {
            auto edges = mesh.getsupersimplices( 1, 0, index );
            for( int e : edges )
                diameter = maximum( diameter, mesh.getDiameter(1,e) );
        }
        assert( diameter > 0 );

        int index_of_entry = volume * subcells_per_simplex + cell;
            
        SparseMatrix::MatrixEntry entry;
        entry.row    = index;
        entry.column = index;
        entry.value  = measure * power_numerical( diameter, -k-ss );
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}

SparseMatrix FEECCechDiffMatrix( const Mesh& mesh, int n, int k )
{
    
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k < n );
    
    // Auxiliary calculations and preparations
    
    const int num_faces = mesh.count_simplices( k   );
    
    const int num_cells = mesh.count_simplices( k+1 );
    
    const int faces_per_cell = k+2; 
    
    const int dim_in      = num_faces;
    const int dim_out     = num_cells;
    const int num_entries = num_cells * faces_per_cell;

    auto face_inclusions = generateSigmas( IndexRange(0,k), IndexRange(0,k+1) );
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_cells;      s++ )
    for( int i = 0; i < faces_per_cell; i++ )
    {
        
        const int face_index = mesh.get_subsimplex( k+1, k, s, i );

        auto std_vector_cellvertices = mesh.getsubsimplices( k+1, 0, s ).getvalues();
        std::sort( std_vector_cellvertices.begin(), std_vector_cellvertices.end() );
        
        auto std_vector_facevertices = mesh.getsubsimplices( k, 0, face_index ).getvalues();
        std::sort( std_vector_facevertices.begin(), std_vector_facevertices.end() );

        int gap_index = 0;
        while( gap_index <= k && std_vector_facevertices[gap_index] == std_vector_cellvertices[gap_index] ) gap_index++;
        
        for( int i = gap_index; i <= k; i++ ) assert( std_vector_facevertices[i] == std_vector_cellvertices[i+1] );
        
        int index_of_entry = s * faces_per_cell + i;
            
        SparseMatrix::MatrixEntry entry;
        entry.row    = s;
        entry.column = face_index;
        entry.value  = sign_power( gap_index );
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}




