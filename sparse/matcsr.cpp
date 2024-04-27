
#include <cmath>
#include <ostream>
#include <utility>
#include <vector>

#include <set>
#include <numeric>

#include "matcsr.hpp"

const bool csr_matrix_verbosity = false;

MatrixCSR::MatrixCSR( 
    int rows,
    int columns,
    const std::vector<int>& A, 
    const std::vector<int>& C, 
    const std::vector<Float>& V
): LinearOperator( rows, columns ),
   A(A), C(C), V(V) 
{
    MatrixCSR::check();
}


MatrixCSR::MatrixCSR( 
    int rows,
    int columns,
    const std::vector<int>&& A, 
    const std::vector<int>&& C, 
    const std::vector<Float>&& V
): LinearOperator( rows, columns ),
   A(std::move(A)), C(std::move(C)), V(std::move(V)) 
{
    MatrixCSR::check();
}


MatrixCSR::MatrixCSR( 
    const SparseMatrix& matrix
): LinearOperator( matrix.getdimout(), matrix.getdimin() ),
   A(0), C(0), V(0) 
{
    matrix.check();
    
    if( csr_matrix_verbosity ) LOG << "Sorting CCO -> CSR: " << matrix.getnumberofentries() << nl;

    if( not matrix.is_sorted() ) {
        matrix.sortandcompressentries();
    }
    
    if( csr_matrix_verbosity ) LOG << "Allocating CCO -> CSR: " << matrix.getnumberofentries() << nl;
    
    const int rows       = matrix.getdimout();
    const int columns    = matrix.getdimin();
    const int numentries = matrix.getnumberofentries();

    // std::vector<int>   A( rows+1,     0  );
    // std::vector<int>   C( numentries, 0  );
    // std::vector<Float> V( numentries, 0. );
    
    A.resize( rows+1     );
    C.resize( numentries );
    V.resize( numentries );
    
    for( int i = 0; i < matrix.getnumberofentries(); i++ ){
        A[ matrix.getentry(i).row+1 ] += 1;
        C[i] = matrix.getentry(i).column;
        V[i] = matrix.getentry(i).value;
    }
 
    for( int i = 1; i < A.size(); i++ ){
        A[i] += A[i-1];
    }

    if( csr_matrix_verbosity ) LOG << "DONE CCO -> CSR: " << matrix.getnumberofentries() << nl;

    for( int i = 0; i < matrix.getnumberofentries(); i++ ){
        assert( 0 <= C[i] and C[i] < columns );
    }
 

    MatrixCSR::check();
    
}

MatrixCSR::MatrixCSR( int rows, int columns )
: MatrixCSR( SparseMatrix( rows, columns ) )
{
    MatrixCSR::check();
}




MatrixCSR::~MatrixCSR()
{
    MatrixCSR::check();
}





MatrixCSR::MatrixCSR( const MatrixCSR& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ),
  A( mat.A ),
  C( mat.C ),
  V( mat.V )
{
    LOG << "*************************************************\n";
    LOG << "*********** WARNING: DEEP COPY ******************\n";
    LOG << "***********  OF CSR MATRIX     ******************\n";
    LOG << "*************************************************\n";
    MatrixCSR::check();
}

MatrixCSR& MatrixCSR::operator=( const MatrixCSR& mat )
{
    LOG << "*************************************************\n";
    LOG << "********** WARNING: DEEP ASSIGN *****************\n";
    LOG << "**********    OF CSR MATRIX     *****************\n";
    LOG << "*************************************************\n";
    assert( getdimin() == mat.getdimin() );
    assert( getdimout() == mat.getdimout() );
    this->A = mat.A;
    this->C = mat.C;
    this->V = mat.V;
    check();
    return *this;
}

MatrixCSR::MatrixCSR( MatrixCSR&& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ),
  A( std::move(mat.A) ),
  C( std::move(mat.C) ),
  V( std::move(mat.V) )
{
    MatrixCSR::check();
}

MatrixCSR& MatrixCSR::operator=( MatrixCSR&& mat )
{
    assert( getdimin() == mat.getdimin() );
    assert( getdimout() == mat.getdimout() );
    this->A = std::move( mat.A );
    this->C = std::move( mat.C );
    this->V = std::move( mat.V );
    check();
    return *this;
}
















void MatrixCSR::check() const
{
    LinearOperator::check();
    
    // check array dimensions and some fix values
    assert( A.size() == getdimout()+1 );
    assert( C.size() == V.size() );
    
    // check that A is ascending 
    for( int p = 1; p <= getdimout(); p++ ) assert( A[p-1] <= A[p] );
    
    // check the final values of A and the validity of the values in C and V
    assert( A[ getdimout() ] == V.size() );
    assert( A[ getdimout() ] == C.size() );
    for( int i = 0; i < C.size(); i++ ) assert( 0 <= C[i] && C[i] < getdimin() && std::isfinite( V[i] ) );

}

std::string MatrixCSR::text() const
{
    std::string str_A, str_C, str_V; 

    for( int i = 0; i < A.size(); i++ ) str_A += ( std::to_string(A[i]) + " " );
    for( int i = 0; i < C.size(); i++ ) str_C += ( std::to_string(C[i]) + " " );
    for( int i = 0; i < V.size(); i++ ) str_V += ( std::to_string(V[i]) + " " );
    
    return std::string("CSRMatrix ") + std::to_string(getdimout()) + "x" + std::to_string(getdimin())
                        + "\nA: " + str_A
                        + "\nC: " + str_C
                        + "\nV: " + str_V;
                        // + "\n";
}

// void MatrixCSR::printplain( std::ostream& os ) const
// {
//     print( os );
// }

void MatrixCSR::apply( FloatVector& dest, const FloatVector& add, Float scaling ) const
{
    check();
    add.check();
    dest.check();
    
    assert( getdimin() == add.getdimension() );
    assert( getdimout() == dest.getdimension() );
    assert( &dest != &add );

    dest.zero();

    Float*       __restrict p_dest = dest.raw();
    const Float* __restrict p_add  = add.raw();
    const int*   __restrict p_A    = A.data();
    const int*   __restrict p_C    = C.data();
    const Float* __restrict p_V    = V.data();
    const int N = A.size();

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int i = 0; i < N-1; i++ ){
        for( int j = p_A[i]; j < p_A[i+1] ; j++ ){
            p_dest[i] += scaling * p_V[j] * p_add[ p_C[j] ];
        }
    }

}








void MatrixCSR::scale ( Float s )
{
    for( auto& v : this->V ) v *= s;
}

bool MatrixCSR::isfinite() const 
{
    for( const Float& value : V )
        if( not std::isfinite(value) ) 
            return false;
    return true;
}

FloatVector MatrixCSR::getDiagonal() const
{ 
    check();
    assert( getdimin() == getdimout() );
    auto ret = FloatVector( getdimin(), 0. );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < getdimout(); r++ ) {
        
        ret[r] = 0.;
        
        for( int i = A[r]; i < A[r+1]; i++ )
            if( C[ i ] == r )
                ret[r] += V[ i ];
        
    }
    
    return ret;
}


static void sort_and_compress_csrdata( std::vector<int>& A, std::vector<int>& C, std::vector<Float>& V );

MatrixCSR MatrixCSR::getTranspose() const 
{
    // gather relevant data
    const int mat_rows = getdimout();
    const int mat_cols = getdimin();
    
    const int*   matA = getA();
    const int*   matC = getC();
    const Float* matV = getV();
    
    const int num_entries = getnumberofentries();
    
    std::vector<int>   B( mat_cols + 1, 0  );
    std::vector<int>   D( num_entries,  0  );
    std::vector<Float> W( num_entries,  0. );
    
    std::vector<int>   Z( mat_cols, 0  );
    
    for( int i = 0; i < num_entries; i++ )
        Z[ matC[i] ]++;

    for( int c = 1; c <= mat_cols; c++ )
        B[c] = B[c-1] + Z[c-1];

    assert( B[mat_cols] == num_entries );
    for( int c = 1; c <= mat_cols; c++ ) assert( B[c-1] <= B[c] );
    

    for( int r = 0; r < mat_rows; r++ ) {
        
        for( int i = matA[r]; i < matA[r+1]; i++ ) {
            
            int   c = matC[i];
            Float v = matV[i];

            int base = B[c];
            int index = --Z[c];
            D[ base + index ] = r;
            W[ base + index ] = v;

            assert( index  >= 0 );
            assert( index  <= B[c+1] - B[c] );
            assert( B[c]   <= base + index );
            assert( B[c+1] >= base + index );
            
        }

    }

    assert( B[mat_cols] == num_entries );
    assert( B[0] == 0 );
    for( int c = 1; c <= mat_cols; c++ ) assert( B[c-1] <= B[c] );
    
    
    // computations done, create matrix 

    sort_and_compress_csrdata(B, D, W );

    return MatrixCSR( mat_cols, mat_rows, std::move(B), std::move(D), std::move(W) );
}










const int*   MatrixCSR::getA() const 
{ 
    return A.data();
}

const int*   MatrixCSR::getC() const 
{ 
    return C.data();
}

const Float* MatrixCSR::getV() const
{ 
    return V.data();
}

        
int MatrixCSR::getnumberofentries() const 
{
    return SIZECAST( V.size() );
}

int MatrixCSR::getnumberofzeroentries() const
{ 
    check();
    
    int ret = 0;
    
    for( const auto& value : V )
        if( value == 0. )
            ret++;
    
    return ret;

}

int MatrixCSR::getmaxrowwidth() const
{
    check();    

    int ret = 0;

    for( int r = 0; r < getdimout(); r++ ) {
        int width = A[r+1] - A[r];
        if( width > ret ) ret = width;
    }

    return ret;
}




// void MatrixCSR::sortentries() const
// {
//     check();
// 
//     // for each row, perform a quick sort of the columns
//     // the following uses an inefficient version of selection sort
//     for( int r = 0; r < A.size()-1; r++ )
//         for( int c1 = A[r]; c1 < A[r+1]; c1++ )
//         for( int c2 = c1+1; c2 < A[r+1]; c2++ )
//             if( C[c1] < C[c2] ) {
//                 std::swap( C[c1], C[c2] );
//                 std::swap( V[c1], V[c2] );
//             }
//     check();
// }

Float MatrixCSR::eigenvalueupperbound() const 
{
    Float ret = 0.;
    for( int r = 0; r < A.size()-1; r++ ) {
        Float candidate_ret = 0.;
        for( int c = A[r]; c < A[r+1]; c++ )
            candidate_ret += absolute( V[c] );
        ret = maximum( ret, candidate_ret );
    }
    return ret;
}



/* Memory size */
        
std::size_t MatrixCSR::memorysize() const
{
    std::size_t ret = 0;
    ret += sizeof(*this);
    ret += A.size() * sizeof(A[0]);
    ret += C.size() * sizeof(C[0]);
    ret += V.size() * sizeof(V[0]);
    return ret;
}







UNUSED static void sort_and_compress_csrdata_perform( std::vector<int>& A, std::vector<int>& C, std::vector<Float>& V )
{

    // return ; 

    const int num_rows = A.size()-1;
    
    assert( A.back() == C.size() );
    assert( A.back() == V.size() );
    
    std::vector<int> nnz( num_rows );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < num_rows; r++ ) {
        
        for( int i = A[r]; i < A[r+1]; i++ ) 
        for( int j = i+1; j < A[r+1]; j++ ) 
        {
            
            // if the columns are the same, first merge the entries
            // there is nothing more to be done
            if( C[i] == C[j] ) {
                V[i] += V[j];
                V[j] = 0.;
                continue;
            }
            
            // if the columns are in the wrong order, 
            // then swap 
            if( C[i] > C[j] ) {
                std::swap( C[i], C[j] );
                std::swap( V[i], V[j] );
            }
            
        }
        
        // swap all the zeroes to the very end
        for( int i = A[r]  ; i < A[r+1]; i++ ) // bubble sort
        for( int j = A[r]+1; j < A[r+1]; j++ ) 
        {
            if( V[j-1] == 0. and V[j] != 0. ) {
                std::swap( C[j-1], C[j] );
                std::swap( V[j-1], V[j] );
            }

        }
        
        // for( int i = A[r]; i < A[r+1]; i++ ) 
        //     LOG << i << space << C[i] << space << V[i] << nl;
            
        for( int i = A[r]+1; i < A[r+1]; i++ ) 
            Assert( C[i-1] < C[i] or V[i] == 0., i, C[i-1], C[i], V[i] );
        
        // Find the first zero within the current row data 
        int first_zero = A[r];
        for(; first_zero < A[r+1] and V[first_zero] != 0.; first_zero++ ) 
        
        for( int i = A[r]; i < first_zero-1; i++ ) Assert( C[i] < C[i+1], i, C[i], C[i+1], V[i], V[i+1] );        
        for( int i = A[r]; i < first_zero;   i++ ) Assert( V[i] != 0., i, C[i], V[i] );
        for( int i = first_zero; i < A[r+1]; i++ ) Assert( V[i] == 0., i, C[i], V[i] );

        nnz[r] = first_zero - A[r]; // we save the number of non-zeroes

        Assert( nnz[r] <= A[r+1] - A[r] );

    }

    std::vector<int> newA( num_rows + 1 );
    
    newA[0] = 0;
    for( int r = 0; r < num_rows; r++ ){
        newA[r+1] = newA[r] + nnz[r];
    }

    for( int r = 0; r < num_rows; r++ ) {
        Assert( nnz[r] == newA[r+1] - newA[r] );
        Assert( newA[r+1] >= newA[r] );
    }
        

    std::vector<int>   newC( newA[ num_rows ] );
    std::vector<Float> newV( newA[ num_rows ] );
    
    // Fill in the new data 
    // #if defined(_OPENMP)
    // #pragma omp parallel for
    // #endif
    for( int r = 0; r < num_rows; r++ ) {

        Assert( nnz[r] <= A[r+1] - A[r] );
        Assert( nnz[r] == newA[r+1] - newA[r] );

        for( int i =          A[r]; i < A[r] + nnz[r]; i++ ) Assert( V[i] != 0. );
        for( int i = A[r] + nnz[r]; i <        A[r+1]; i++ ) Assert( V[i] == 0. );
        
        for( int i = 0; i < nnz[r]; i++ ) {
            newC[ newA[r] + i ] = C[ A[r] + i ];
            newV[ newA[r] + i ] = V[ A[r] + i ];
        }

    }

    nnz.clear(); // NOTE: just to free some memory


    for( int r = 0; r < num_rows; r++ )
    for( int i = newA[r]+1; i < newA[r+1]; i++ ) 
        assert( newC[i-1] < newC[i]);
    
    for( int i = 0; i < newV.size(); i++ ) 
        assert( newV[i] != 0. );

    A = std::move( newA );
    C = std::move( newC );
    V = std::move( newV );
    
}





static void sort_and_compress_csrdata_reduced( std::vector<int>& A, std::vector<int>& C, std::vector<Float>& V )
{

    // return ; 

    const int num_rows = A.size()-1;
    
    assert( A.back() == C.size() );
    assert( A.back() == V.size() );
    
    // sort the row data and count the zero entries
    std::vector<int> nnz( num_rows );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < num_rows; r++ ) {
        
        for( int i = A[r]; i < A[r+1]; i++ ) 
        for( int j = i+1; j < A[r+1]; j++ ) 
        {
            
            // if the columns are the same, first merge the entries
            // there is nothing more to be done
            if( C[i] == C[j] ) {
                V[i] += V[j];
                V[j] = 0.;
                continue;
            }
            
            // if the columns are in the wrong order, 
            // then swap 
            if( C[i] > C[j] ) {
                std::swap( C[i], C[j] );
                std::swap( V[i], V[j] );
            }
            
        }
        
        // swap all the zeroes to the very end
        for( int i = A[r]  ; i < A[r+1]; i++ ) // bubble sort
        for( int j = A[r]+1; j < A[r+1]; j++ ) 
        {
            if( V[j-1] == 0. and V[j] != 0. ) {
                std::swap( C[j-1], C[j] );
                std::swap( V[j-1], V[j] );
            }

        }
        
        // for( int i = A[r]; i < A[r+1]; i++ ) 
        //     LOG << i << space << C[i] << space << V[i] << nl;
            
        for( int i = A[r]+1; i < A[r+1]; i++ ) 
            Assert( C[i-1] < C[i] or V[i] == 0., i, C[i-1], C[i], V[i] );
        
        // Find the first zero within the current row data 
        int first_zero = A[r];
        for(; first_zero < A[r+1] and V[first_zero] != 0.; first_zero++ ) 
        
        for( int i = A[r]; i < first_zero-1; i++ ) Assert( C[i] < C[i+1], i, C[i], C[i+1], V[i], V[i+1] );        
        for( int i = A[r]; i < first_zero;   i++ ) Assert( V[i] != 0., i, C[i], V[i] );
        for( int i = first_zero; i < A[r+1]; i++ ) Assert( V[i] == 0., i, C[i], V[i] );

        nnz[r] = first_zero - A[r]; // we save the number of non-zeroes

        Assert( nnz[r] <= A[r+1] - A[r] );

    }

    
    // Update C and V (sequential)
    int index_target = 0;
    
    for( int r = 0; r < num_rows; r++ )
    {
        for( int i = A[r]; i < A[r] + nnz[r]; i++ )
        {
            assert( index_target <= i );
            C[index_target] = C[i];
            V[index_target] = V[i];
            index_target++;
        }
    }


    // Update A (sequential)
    for( int r = 0; r < num_rows; r++ )
    {
        A[r+1] = nnz[r];
    }

    A[0] = 0;
    for( int r = 1; r <= num_rows; r++ )
    {
        A[r] += A[r-1];
    }

    C.resize( A[num_rows] );
    V.resize( A[num_rows] );

    assert( A[num_rows] == std::accumulate( nnz.begin(), nnz.end(), 0 ) );
}


static void sort_and_compress_csrdata( std::vector<int>& A, std::vector<int>& C, std::vector<Float>& V )
{
    sort_and_compress_csrdata_reduced( A, C, V );
}


void MatrixCSR::compressentries() const
{
    auto& p_A = const_cast< std::vector<int>& >(A);
    auto& p_C = const_cast< std::vector<int>& >(C);
    auto& p_V = const_cast< std::vector<Float>& >(V);

    sort_and_compress_csrdata(p_A,p_C,p_V);
}



MatrixCSR MatrixCSRAddition( const MatrixCSR& mat1, const MatrixCSR& mat2, Float s1, Float s2 )
{
    // gather relevant data
    const int mat1_rows = mat1.getdimout();
    const int mat1_cols = mat1.getdimin();
    const int mat2_rows = mat2.getdimout();
    const int mat2_cols = mat2.getdimin();
    
    const int* mat1A = mat1.getA();
    const int* mat2A = mat2.getA();
    
    const int* mat1C = mat1.getC();
    const int* mat2C = mat2.getC();

    const Float* mat1V = mat1.getV();
    const Float* mat2V = mat2.getV();

    Assert( mat1_cols == mat2_cols );
    Assert( mat1_rows == mat2_rows );

    const int matn_rows = mat1_rows;
    const int matn_cols = mat1_cols;

    const int num_entries = mat1A[matn_rows] + mat2A[matn_rows];
    
    std::vector<int>   A( matn_rows + 1, 0 );
    std::vector<int>   C( num_entries, 0  );
    std::vector<Float> V( num_entries, 0. );
    
    A[0] = 0;
    for( int r = 1; r <= matn_rows; r++ )
        A[r] = A[r-1] + (mat1A[r]-mat1A[r-1]) + (mat2A[r]-mat2A[r-1]);
    Assert( A[matn_rows] == num_entries );
    
    for( int r = 0; r < matn_rows; r++ ) {
        
        int i = 0;
        for( int c1 = mat1A[r]; c1 < mat1A[r+1]; c1++ ) C[ A[r] + (i++) ] = mat1C[c1];
        for( int c2 = mat2A[r]; c2 < mat2A[r+1]; c2++ ) C[ A[r] + (i++) ] = mat2C[c2];
        Assert( A[r] + i == A[r+1], r, A[r], i, A[r+1] );

    }
    
    for( int r = 0; r < matn_rows; r++ ) {
        
        int i = 0;
        for( int v1 = mat1A[r]; v1 < mat1A[r+1]; v1++ ) V[ A[r] + (i++) ] = s1 * mat1V[v1];
        for( int v2 = mat2A[r]; v2 < mat2A[r+1]; v2++ ) V[ A[r] + (i++) ] = s2 * mat2V[v2];
        assert( A[r] + i == A[r+1] );

    }
    
    // computations done, create matrix 

    sort_and_compress_csrdata(A, C, V );

    return MatrixCSR( matn_rows, matn_cols, std::move(A), std::move(C), std::move(V) );

}

MatrixCSR MatrixCSRMultiplication( const MatrixCSR& mat1, const MatrixCSR& mat2 )
{
    // gather relevant data
    const int mat1_rows = mat1.getdimout();
    const int mat1_cols = mat1.getdimin();
    const int mat2_rows = mat2.getdimout();
    const int mat2_cols = mat2.getdimin();
    
    const int* mat1A = mat1.getA();
    const int* mat2A = mat2.getA();
    
    const int* mat1C = mat1.getC();
    const int* mat2C = mat2.getC();

    const Float* mat1V = mat1.getV();
    const Float* mat2V = mat2.getV();

    Assert( mat1_cols == mat2_rows );

    const int matn_rows = mat1_rows;
    const int matn_cols = mat2_cols;

    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication: create index list A, allocate C and V" << nl; 
    
    std::vector<int> A( mat1_rows + 1, 0 );

    LOG << mat1.text() << nl << mat2.text() << nl;

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < matn_rows; r++ )
        for( int c1 = mat1A[r];           c1 < mat1A[r+1];           c1++ )
        for( int c2 = mat2A[ mat1C[c1] ]; c2 < mat2A[ mat1C[c1]+1 ]; c2++ )
            A[r+1]++;

    // for( int r = 0; r <= matn_rows; r++ ) LOG << A[r] << space; LOG << nl;
    
    for( int r = 1; r <= matn_rows; r++ ) A[r] += A[r-1];

    // for( int r = 0; r <= matn_rows; r++ ) LOG << A[r] << space; LOG << nl;
    
    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication: Temporary number of elements:" << A[matn_rows] << nl;

    std::vector<int>   C( A[matn_rows] );
    std::vector<Float> V( A[matn_rows] );

    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication: allocated" << nl; 
    
    // compute entries     

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < matn_rows; r++ ) {

        int c_curr = A[r];
        
        for( int c1 = mat1A[r];           c1 < mat1A[r+1];           c1++ )
        for( int c2 = mat2A[ mat1C[c1] ]; c2 < mat2A[ mat1C[c1]+1 ]; c2++ )
        {
            C[ c_curr ] = mat2C[ c2 ];
            V[ c_curr ] = mat1V[ c1 ] * mat2V[ c2 ];
            c_curr++;
        }
            
        Assert( c_curr == A[r+1], r, c_curr, A[r+1] ); 

    }
        
    // computations done
    
    // create matrix 
    
    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication: sort and compress" << nl; 
    
    sort_and_compress_csrdata( A, C, V );

    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication: return" << nl; 
    
    return MatrixCSR( matn_rows, matn_cols, std::move(A), std::move(C), std::move(V) );

}







MatrixCSR MatrixCSRMultiplication_reduced( const MatrixCSR& mat1, const MatrixCSR& mat2 )
{
    // gather relevant data
    const int mat1_rows = mat1.getdimout();
    const int mat1_cols = mat1.getdimin();
    const int mat2_rows = mat2.getdimout();
    const int mat2_cols = mat2.getdimin();
    
    const int* mat1A = mat1.getA();
    const int* mat2A = mat2.getA();
    
    const int* mat1C = mat1.getC();
    const int* mat2C = mat2.getC();

    const Float* mat1V = mat1.getV();
    const Float* mat2V = mat2.getV();

    Assert( mat1_cols == mat2_rows );

    const int matn_rows = mat1_rows;
    const int matn_cols = mat2_cols;

    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication reduced: create index list A and fill in" << nl; 
    
    std::vector<int> A( mat1_rows + 1, 0 );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < matn_rows; r++ )
    {
        std::set<int> column_indices;

        for( int c1 = mat1A[r];           c1 < mat1A[r+1];           c1++ )
        for( int c2 = mat2A[ mat1C[c1] ]; c2 < mat2A[ mat1C[c1]+1 ]; c2++ )
            if( mat1V[ c1 ] != 0. && mat2V[ c2 ] != 0. )
                column_indices.insert( mat2C[c2] );

        A[r+1] = column_indices.size();
    }

    for( int r = 1; r <= matn_rows; r++ ) A[r] += A[r-1];
    for( int r = 1; r <= matn_rows; r++ ) assert( A[r-1] <= A[r] );

    // for( int r = 0; r <= matn_rows; r++ ) LOG << A[r] << space; LOG << nl;    
    // LOG << mat1.text() << nl << mat2.text() << nl;

    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication reduced: Temporary number of elements:" << A[matn_rows] << nl;

    std::vector<int>   C( A[matn_rows], 0  );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < matn_rows; r++ )
    {
        std::set<int> column_indices;

        for( int c1 = mat1A[r];           c1 < mat1A[r+1];           c1++ )
        for( int c2 = mat2A[ mat1C[c1] ]; c2 < mat2A[ mat1C[c1]+1 ]; c2++ )
            if( mat1V[ c1 ] != 0. && mat2V[ c2 ] != 0. )
                column_indices.insert( mat2C[c2] );

        Assert( A[r] + column_indices.size() == A[r+1] );

        int i = 0;
        for( int c : column_indices ) {
            C[ A[r] + i ] = c;
            i++;
        }
        assert( i == column_indices.size() );

        for( int c = A[r]+1; c < A[r+1]; c++ ){
            assert( C[c] > C[c-1] );
        }
        

        Assert( A[r]+i == A[r+1], A[r], i, A[r+1] );

    }

    std::vector<Float> V( A[matn_rows], 0. );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int r = 0; r < matn_rows; r++ )
    {
        for( int c = A[r]; c < A[r+1]; c++ )
        for( int c1 = mat1A[r];           c1 < mat1A[r+1];           c1++ )
        for( int c2 = mat2A[ mat1C[c1] ]; c2 < mat2A[ mat1C[c1]+1 ]; c2++ )
        {
            assert( A[r] < A[r+1] );
            
            if( C[c] == mat2C[c2] )
                if( mat1V[ c1 ] != 0. && mat2V[ c2 ] != 0. )
                    V[ c ] += mat1V[ c1 ] * mat2V[ c2 ];
        }
        
    }

    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication reduced: compressing" << nl; 
    sort_and_compress_csrdata( A, C, V );
    
    assert( A.back() == A[matn_rows] );
    if( csr_matrix_verbosity ) LOG << "MatrixCSRMultiplication reduced: return " << A.back() << nl; 
    
    return MatrixCSR( matn_rows, matn_cols, std::move(A), std::move(C), std::move(V) );

}







DiagonalOperator InverseDiagonalPreconditioner( const MatrixCSR& mat )
{ 
    mat.check();
    assert( mat.getdimin() == mat.getdimout() );

    auto diag = mat.getDiagonal();
    
    for( int r = 0; r < mat.getdimin(); r++ )
    {
        assert( diag[r] >= 0. );
        
        if( diag[r] > 0. ) diag[r] = 1. / diag[r];
    }
    
    return DiagonalOperator( diag );
    
//     auto ret = FloatVector( mat.getdimin(), 0. );
// 
//    #if defined(_OPENMP)
//    #pragma omp parallel for
//    #endif
//    for( int r = 0; r < mat.getdimin(); r++ ) {
//         
//         ret[r] = 0.;
//         
//         for( int c = mat.A[r]; c < mat.A[r+1]; c++ )
//             if( mat.C[ c ] == r )
//                 ret[r] += mat.V[ c ];
//     
//         assert( ret[r] >= 0. );
//         
//         if( ret[r] > 0. ) ret[r] = 1. / ret[r];
//         
//     }
//      
}



































