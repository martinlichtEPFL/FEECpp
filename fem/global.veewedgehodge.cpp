
#include <algorithm>
#include <vector>
#include <tuple>

#include "../basic.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/global.veewedgehodge.hpp"



SparseMatrix FEECBrokenVeeMatrix( const Mesh& mesh, int n, int k, int r , int l, int s, FloatVector field )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( r >= 0 );
    assert( s >= 0 );
    assert( k >= 0 && k <= n );
    assert( l >= 0 && l <= n );
    assert( k >= l );
    assert( field.isfinite() );

    assert( k == l ); // restricted special case for now...
    
    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim_in  = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int localdim_out = binomial_integer( n+(r+s), n ) * binomial_integer( n+1, (k-l) );
    
    const auto mis_input  = generateMultiIndices( IndexRange(0,n), r   );
    const auto mis_factor = generateMultiIndices( IndexRange(0,n), s   );
    const auto mis_output = generateMultiIndices( IndexRange(0,n), r+s );

    std::vector<std::array<int,3>> couplings;

    for( int m_i = 0; m_i < mis_input.size();  m_i++ )
    for( int m_f = 0; m_f < mis_factor.size(); m_f++ )
    {
        const auto& mi_input  = mis_input[m_i];
        const auto& mi_factor = mis_factor[m_f];
        
        const auto mi_prod = mi_input + mi_factor;
        
        for( int m_o = 0; m_o < mis_output.size(); m_o++ )
        {
            if( mi_prod == mis_output[m_o] ) {
                couplings.push_back( {m_i,m_f,m_o} );
            }                
        }
    }

    // LOG << "Couplings: " << couplings.size() << nl; for( const auto& a : couplings ) LOG << a[0] << space << a[1] << space << a[2] << nl;


    
    const int dim_in      = num_simplices * localdim_in;
    const int dim_out     = num_simplices * localdim_out;
    const int num_entries = num_simplices * couplings.size() * binomial_integer( n+1, k );

    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        DenseMatrix GPM    = mesh.getGradientProductMatrix( n, s );
        DenseMatrix formMM = SubdeterminantMatrix( GPM, k );
        assert( formMM.isfinite() );

        assert( not formMM.iszero() );
        assert( formMM.issquare() );
        assert( formMM.issymmetric() );

        for( int coupling_index = 0; coupling_index < couplings.size();        coupling_index++ ) 
        for( int f_i            = 0;            f_i < binomial_integer(n+1,k); f_i++            )
        {

            int m_i = couplings[coupling_index][0];
            int m_f = couplings[coupling_index][1];
            int m_o = couplings[coupling_index][2];

            SparseMatrix::MatrixEntry entry;
            entry.row    = s * mis_output.size() + m_o;
            entry.column = s * mis_input.size() * binomial_integer(n+1,k)  + m_i * binomial_integer(n+1,k) + f_i;
            
            entry.value = 0.;
            for( int g_i = 0; g_i < binomial_integer(n+1,k); g_i++ )
                entry.value += field[ s * mis_factor.size() * binomial_integer(n+1,k)  + m_f * binomial_integer(n+1,k) + g_i ] * formMM( g_i, f_i );

            int index_of_entry = s * couplings.size() * binomial_integer( n+1, k ) + coupling_index * binomial_integer( n+1, k ) + f_i;
            
            ret.setentry( index_of_entry,  entry );
            assert( std::isfinite(entry.value) );
        }        
        
    }
    
    assert( ret.isfinite() );

    return ret;
}













SparseMatrix FEECBrokenWedgeMatrix( const Mesh& mesh, int n, int k, int r , int l, int s, FloatVector field )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( r >= 0 );
    assert( s >= 0 );
    assert( k >= 0 && k <= n );
    assert( l >= 0 && l <= n );
    assert( field.isfinite() );

    // Auxiliary calculations and preparations
    
    const int num_simplices = mesh.count_simplices( n );
        
    const int localdim_in  = binomial_integer( n+r, n ) * binomial_integer( n+1, k );
    
    const int localdim_out = binomial_integer( n+(r+s), n ) * binomial_integer( n+1, k+l );
    
    const auto mis_input  = generateMultiIndices( IndexRange(0,n), r   );
    const auto mis_factor = generateMultiIndices( IndexRange(0,n), s   );
    const auto mis_output = generateMultiIndices( IndexRange(0,n), r+s );

    std::vector<std::array<int,3>> mi_couplings;

    for( int m_i = 0; m_i < mis_input.size();  m_i++ )
    for( int m_f = 0; m_f < mis_factor.size(); m_f++ )
    {
        const auto& mi_input  = mis_input[m_i];
        const auto& mi_factor = mis_factor[m_f];
        
        const auto mi_prod = mi_input + mi_factor;
        
        for( int m_o = 0; m_o < mis_output.size(); m_o++ )
        {
            if( mi_prod == mis_output[m_o] ) {
                mi_couplings.push_back( {m_i,m_f,m_o} );
            }                
        }
    }


    const auto sigmas_input  = generateSigmas( IndexRange(1,k  ), IndexRange(0,n) );
    const auto sigmas_factor = generateSigmas( IndexRange(1,  l), IndexRange(0,n) );
    const auto sigmas_output = generateSigmas( IndexRange(1,k+l), IndexRange(0,n) );

    std::vector<std::array<int,4>> sigma_couplings;

    for( int s_i = 0; s_i < sigmas_input.size();  s_i++ )
    for( int s_f = 0; s_f < sigmas_factor.size(); s_f++ )
    {
        const auto& sigma_input  = sigmas_input[s_i];
        const auto& sigma_factor = sigmas_factor[s_f];
        
        int signum;
        IndexMap sigma_prod = mergeSigmas( sigma_factor, sigma_input, signum );

        if( signum == 0 ) continue;
        
        for( int s_o = 0; s_o < sigmas_output.size(); s_o++ )
        {
            if( sigma_prod == sigmas_output[s_o] ) {
                sigma_couplings.push_back( {s_i,s_f,s_o,signum} );
            }                
        }
    }

    // LOG << "Couplings: " << couplings.size() << nl; for( const auto& a : couplings ) LOG << a[0] << space << a[1] << space << a[2] << nl;


    
    const int dim_in      = num_simplices * localdim_in;
    const int dim_out     = num_simplices * localdim_out;
    const int num_entries = num_simplices * mi_couplings.size() * sigma_couplings.size();

    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        int local_index_of_entry = 0;

        for( int mi_coupling_index = 0;    mi_coupling_index < mi_couplings.size();       mi_coupling_index++    ) 
        for( int sigma_coupling_index = 0; sigma_coupling_index < sigma_couplings.size(); sigma_coupling_index++ ) 
        {

            int m_i = mi_couplings[mi_coupling_index][0];
            int m_f = mi_couplings[mi_coupling_index][1];
            int m_o = mi_couplings[mi_coupling_index][2];

            int s_i = sigma_couplings[sigma_coupling_index][0];
            int s_f = sigma_couplings[sigma_coupling_index][1];
            int s_o = sigma_couplings[sigma_coupling_index][2];

            int signum = sigma_couplings[sigma_coupling_index][3];

            SparseMatrix::MatrixEntry entry;
            entry.row    =                s * mis_output.size() * binomial_integer(n+1,k+l) + m_o * binomial_integer(n+1,k+l) + s_o;
            entry.column =                s * mis_input.size()  * binomial_integer(n+1,k  ) + m_i * binomial_integer(n+1,k  ) + s_i;
            entry.value = signum * field[ s * mis_factor.size() * binomial_integer(n+1,  l) + m_f * binomial_integer(n+1,  l) + s_f ];
            
            assert( std::isfinite(entry.value) );

            int index_of_entry = s * mi_couplings.size() * sigma_couplings.size() + local_index_of_entry;
            
            ret.setentry( index_of_entry, entry );
            local_index_of_entry++;
        }        
        
        assert( local_index_of_entry == mi_couplings.size() * sigma_couplings.size() );
    }
    
    assert( ret.isfinite() );

    return ret;
}














SparseMatrix FEECBrokenHodgeMatrix( const Mesh& mesh, int n, int k, int r )
{
    
    assert( r >= 0 );
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( k >= 0 && k <= n );
    
    
    const int num_simplices = mesh.count_simplices( n );

    auto sigmas_input  = generateSigmas( IndexRange(1,k  ), IndexRange(0,n) );
    auto sigmas_test   = generateSigmas( IndexRange(1,k  ), IndexRange(0,n) );
    auto applejack_test = std::remove_if( sigmas_test.begin(), sigmas_test.end(),
                    [k](IndexMap im) -> bool { assert( im.isstrictlyascending() ); return   k != 0 and im[1] == 0; }
    );
    sigmas_test.erase( applejack_test, sigmas_test.end() );

    auto sigmas_output = generateSigmas( IndexRange(1,n-k), IndexRange(0,n) );
    auto sigmas_hodge  = generateSigmas( IndexRange(1,n-k), IndexRange(0,n) );
    auto applejack_hodge = std::remove_if( sigmas_hodge.begin(), sigmas_hodge.end(),
                    [n,k](IndexMap im) -> bool { assert( im.isstrictlyascending() ); return n-k != 0 and im[1] == 0; }
    );
    sigmas_hodge.erase( applejack_hodge, sigmas_hodge.end() );

    // for( auto x : sigmas_hodge ) LOG << x << nl;
    
    const int dim_polynomials = binomial_integer( n+r, n   );
    const int dim_sigmas_in   = binomial_integer( n+1,   k );
    const int dim_sigmas_out  = binomial_integer( n+1, n-k );
    
    assert( sigmas_input.size()  == dim_sigmas_in       );
    assert( sigmas_output.size() == dim_sigmas_out      );
    Assert( sigmas_test.size()   == sigmas_hodge.size(), n, k, sigmas_test.size(), sigmas_hodge.size() );
    
    
    // 1. Assemble the algebraic matrix

    DenseMatrix wedge_matrix( sigmas_test.size(), sigmas_hodge.size(), 0. );

    const IndexMap volume_form( IndexRange(1,n), IndexRange(0,n), [](int i)->int{return i;});

    for( int s_h = 0; s_h < sigmas_hodge.size(); s_h++ )
    for( int s_t = 0; s_t < sigmas_test.size();  s_t++ )
    {
        const auto& sigma_hodge = sigmas_hodge[s_h];
        const auto& sigma_test  = sigmas_test[s_t];
        
        int signum;
        IndexMap sigma_prod = mergeSigmas( sigma_test, sigma_hodge, signum );

        if( signum == 0 ) continue;

        // LOG << sigma_prod << nl << volume_form << nl;

        assert( sigma_prod.isstrictlyascending() );

        assert( sigma_prod == volume_form );

        wedge_matrix( s_t, s_h ) = signum;
    }

    assert( wedge_matrix.issquare() );

    const DenseMatrix wedge_matrix_inv = Inverse( wedge_matrix ); // TODO: This inversion can be made much simpler ... 

    assert( ( wedge_matrix_inv * wedge_matrix - IdentityMatrix( wedge_matrix.getdimin() ) ).is_numerically_small() );

    DenseMatrix wedge_matrix_inv_full( sigmas_output.size(), sigmas_input.size(), 0. );
    
    for( int s_h = 0; s_h < sigmas_hodge.size();  s_h++ )
    for( int s_t = 0; s_t < sigmas_test.size();   s_t++ )
    for( int s_i = 0; s_i < sigmas_input.size();  s_i++ )
    for( int s_o = 0; s_o < sigmas_output.size(); s_o++ )
    {
        if( sigmas_test[s_t] != sigmas_input[s_i] or sigmas_hodge[s_h] != sigmas_output[s_o] ) {
            continue;
        } else {
            wedge_matrix_inv_full( s_o, s_i ) = wedge_matrix_inv( s_h, s_t );
        }
    }

    assert( wedge_matrix_inv_full.isfinite() );
    // for( auto x : sigmas_output ) LOG << x << nl;
    // for( auto x : sigmas_input ) LOG << x << nl;
    // LOG << nl << wedge_matrix_inv_full << nl;
    
    
    
    
    // 2. Run over the simplices and compile the local entries 
    
    // const int localdim_in  = dim_polynomials * dim_sigmas_in;
    // const int localdim_out = dim_polynomials * dim_sigmas_out;
    
    const int dim_in      = num_simplices * dim_polynomials * dim_sigmas_in;
    const int dim_out     = num_simplices * dim_polynomials * dim_sigmas_out;
    const int num_entries = num_simplices * dim_polynomials * dim_sigmas_in * dim_sigmas_out;
    
    SparseMatrix ret( dim_out, dim_in, num_entries );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_simplices; s++ )
    {
        
        const DenseMatrix GPM               = mesh.getGradientProductMatrix( n, s );
        
        const DenseMatrix formMM            = SubdeterminantMatrix( GPM, k );

        assert( formMM.isfinite() );
        
        const DenseMatrix full_local_matrix = wedge_matrix_inv_full * formMM;

        assert( full_local_matrix.isfinite() );
        
        const Float scaling = mesh.getMeasure( n, s ) * factorial_integer(n);
    
        for( int p   = 0; p   < dim_polynomials; p++   )
        for( int s_i = 0; s_i < dim_sigmas_in;   s_i++ )
        for( int s_o = 0; s_o < dim_sigmas_out;  s_o++ )
        {

            int index_of_entry = s * dim_polynomials * dim_sigmas_in * dim_sigmas_out + p * dim_sigmas_in * dim_sigmas_out + s_i * dim_sigmas_out + s_o;
            
            SparseMatrix::MatrixEntry entry;
            entry.row    = s * dim_polynomials * dim_sigmas_out + p * dim_sigmas_out + s_o;
            entry.column = s * dim_polynomials * dim_sigmas_in  + p * dim_sigmas_in  + s_i;
            entry.value  = full_local_matrix( s_o, s_i ) * scaling;
                
            ret.setentry( index_of_entry, entry );

        }
        
    }
    
    return ret;
}
