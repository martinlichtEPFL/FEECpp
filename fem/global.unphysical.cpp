
#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../operators/linearoperator.hpp"
#include "../mesh/mesh.hpp"
#include "../utility/random.hpp"

#include "indexfunctions.hpp"

#include "../fem/global.unphysical.hpp"




SparseMatrix FEECCanonicalizeBroken( const Mesh& mesh, int n, int k, int r )
{
    Assert( 0 <= n );
    Assert( n <= mesh.getinnerdimension() );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );


    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int poly_size = binomial_integer( n+r, n );
    const int form_size = binomial_integer( n+1, k );

    const int dim_in      = num_simplices * poly_size * form_size;
    const int dim_out     = num_simplices * poly_size * form_size;
    const int num_entries = num_simplices * poly_size * form_size * form_size;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    
    // Calculate local matrix 

    DenseMatrix Aux1( n+1, n+1, 0. );
    for( int i = 1; i <= n; i++ ) {
        Aux1(i,i) = 1.;
        Aux1(i,0) = -1.;
    }
    // DenseMatrix Aux1 = IdentityMatrix(n+1) - DenseMatrix( n+1, n+1, 1./(n+1) );
    
    const DenseMatrix Aux = SubdeterminantMatrix( Aux1, k );

    assert( Aux.issquare() and Aux.getdimout() == form_size );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s  = 0;  s < num_simplices;  s++ )
    for( int p  = 0;  p <     poly_size;  p++ )
    for( int f1 = 0; f1 <     form_size; f1++ )
    for( int f2 = 0; f2 <     form_size; f2++ )
    {
        int index_of_entry = s * poly_size * form_size * form_size + p * form_size * form_size + f1 * form_size + f2;

        SparseMatrix::MatrixEntry entry;
        entry.row    = s * poly_size * form_size + p * form_size + f1;
        entry.column = s * poly_size * form_size + p * form_size + f2;
        entry.value  = Aux( f1, f2 );
        
        ret.setentry( index_of_entry, entry );
    }
    
    return ret;

}





SparseMatrix FEECRandomizeBroken( const Mesh& mesh, int n, int k, int r, Float base_alpha )
{
    Assert( 0 <= n );
    Assert( n <= mesh.getinnerdimension() );
    Assert( 0 <= r );
    Assert( 0 <= k and k <= n );


    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int poly_size = binomial_integer( n+r, n );
    const int form_size = binomial_integer( n+1, k );

    const int dim_in      = num_simplices * poly_size * form_size;
    const int dim_out     = num_simplices * poly_size * form_size;
    const int num_entries = num_simplices * poly_size * form_size * form_size;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    
    std::vector<DenseMatrix> auxiliaries( n+1, DenseMatrix(form_size,form_size,notanumber) );

    // Calculate local matrix 
    for( int t = 0; t <= n; t++ )
    {
//         Float alpha = ( 0 <= base_alpha and base_alpha <= 1.0 ) ? base_alpha : random_uniform();
        Float alpha = 1.;
        Assert( 0. <= alpha and alpha <= 1. );

        DenseMatrix Aux1( n+1, n+1, 0. );
        for( int i = 0; i < t; i++ ) {
            Aux1(i,i) =  1.;
            Aux1(i,t) = -alpha;
        }
        Aux1(t,t) = 1. - alpha;    
        for( int i = t+1; i <= n; i++ ) {
            Aux1(i,i) = 1.;
            Aux1(i,t) = -alpha;
        }
        
        const DenseMatrix Aux2 = SubdeterminantMatrix( Aux1, k );

        assert( Aux2.issquare() and Aux2.getdimout() == form_size );        

        auxiliaries[t] = Aux2;
    }
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s  = 0;  s < num_simplices;  s++ )
    {
        int t = random_integer() % (n+1);
        assert( 0 <= t and t <= n );
        
        for( int p  = 0;  p <     poly_size;  p++ )
        for( int f1 = 0; f1 <     form_size; f1++ )
        for( int f2 = 0; f2 <     form_size; f2++ )
        {
            int index_of_entry = s * poly_size * form_size * form_size + p * form_size * form_size + f1 * form_size + f2;
            
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * poly_size * form_size + p * form_size + f1;
            entry.column = s * poly_size * form_size + p * form_size + f2;
            entry.value  = auxiliaries[t]( f1, f2 );
            
            ret.setentry( index_of_entry, entry );
        }
    }
    
    return ret;

}
