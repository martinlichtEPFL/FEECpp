
#include "functions.hpp"

#include <new>
#include <vector>


#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../combinatorics/heappermgen.hpp"

#include "densematrix.hpp"




/************************************
*****                          ******
*****    TRANSPOSE             ******
*****                          ******
************************************/



DenseMatrix Transpose( const DenseMatrix& src ) 
{
    src.check();
    DenseMatrix ret( src.getdimin(), src.getdimout() );
    for( int r = 0; r < src.getdimout(); r++ )
    for( int c = 0; c < src.getdimin();  c++ )
        ret(c,r) = src(r,c);
    ret.check();
    return ret;
}

DenseMatrix TransposeSquare( const DenseMatrix& src ) 
{
    src.check();
    return Transpose( src );
}

// void TransposeInSitu( DenseMatrix& src )
// {
//     src.check();
//     const int numrows = src.getdimout();
//     const int numcols = src.getdimin();
// 
//     for( int start = 0; start < numcols * numrows; ++start )
//     {
// 
//         int next = start;
//         int i = 0;
// 
//         do {
//             ++i;
//             next = (next % numrows) * numcols + next / numrows;
//         } while ( next > start );
// 
//         if( next >= start && i != 1 )
//         {
//             const Float tmp = src( start / numcols, start % numcols );
//             next = start;
//             do {
//                 i = (next % numrows) * numcols + next / numrows;
//                 src( next / numcols, next % numcols ) = ( i == start ) ? tmp : src( i / numcols, i % numcols );
//                 next = i;
//             } while( next > start );
//         }
// 
//     }
// 
//     // We cannot change the dimensions of the matrix `src`
// 
// }

void TransposeSquareInSitu( DenseMatrix& src ) 
{
    src.check();
    for( int r =   0; r < src.getdimout(); r++ )
    for( int c = r+1; c <  src.getdimin(); c++ )
    { 
        Float t = src(c,r); src(c,r) = src(r,c); src(r,c) = t; 
    }
    src.check();
}

//         ret.entries.at( c * src.getdimout() + r ) = src.entries.at( r * src.getdimin() + c );







// DenseMatrix skip_row( int r, const DenseMatrix& mat )
// {
//     assert( 0 <= r and r < mat.getdimout() );
//     DenseMatrix ret( mat.getdimout()-1, mat.getdimin() );
//     for( int i = 0;   i < r;               i++ ) ret.setrow( i,   mat.getrow(i) );
//     for( int i = r+1; i < mat.getdimout(); i++ ) ret.setrow( i-1, mat.getrow(i) );
//     return ret; 
// }

// DenseMatrix skip_column( int c, const DenseMatrix& mat )
// {
//     return Transpose( skip_row(c, Transpose(mat) ) );
// }





/************************************
*****                          ******
*****    DETERMINANT           ******
*****                          ******
************************************/


Float Determinant( const DenseMatrix& A )
{
    assert( A.issquare() );

    // LOG << A.getdimin() << nl;
    
    if( A.getdimin() == 0 ) {
        
        return 1.;
        
    } else if( A.getdimin() == 1 ) {
        
        return A(0,0);
        
    } else if( A.getdimin() == 2 ) {
        
        return A(0,0) * A(1,1) - A(0,1) * A(1,0);
        
    } else if( A.getdimin() == 3 ) {
        
        return + A(0,0) * ( A(1,1)*A(2,2) - A(2,1)*A(1,2) )
               - A(1,0) * ( A(0,1)*A(2,2) - A(2,1)*A(0,2) )
               + A(2,0) * ( A(0,1)*A(1,2) - A(1,1)*A(0,2) )
               ;
        
        // return + A(0,0) * A(1,1) * A(2,2) // 1 2 3 + 
        //        - A(0,0) * A(1,2) * A(2,1) // 1 3 2 - 
        //        - A(0,1) * A(1,0) * A(2,2) // 2 1 3 - 
        //        + A(0,1) * A(1,2) * A(2,0) // 2 3 1 + 
        //        - A(0,2) * A(1,1) * A(2,0) // 3 2 1 - 
        //        + A(0,2) * A(1,0) * A(2,1) // 3 1 2 + 
        //        ;
        
    // } else if( A.getdimin() == 4 ) {
        
    //     return + A(0,0) * A(1,1) * A(2,2) * A(3,3) // 0 1 2 3 + 
    //            - A(0,0) * A(1,2) * A(2,1) * A(3,3) // 0 2 1 3 - 
    //            - A(0,1) * A(1,0) * A(2,2) * A(3,3) // 1 0 2 3 - 
    //            + A(0,1) * A(1,2) * A(2,0) * A(3,3) // 1 2 0 3 + 
    //            - A(0,2) * A(1,1) * A(2,0) * A(3,3) // 2 1 0 3 - 
    //            + A(0,2) * A(1,0) * A(2,1) * A(3,3) // 2 0 1 3 + 
               
    //            - A(0,0) * A(1,1) * A(2,3) * A(3,2) // 0 1 3 2 - 
    //            + A(0,0) * A(1,2) * A(2,3) * A(3,1) // 0 2 3 1 + 
    //            + A(0,1) * A(1,0) * A(2,3) * A(3,2) // 1 0 3 2 + 
    //            - A(0,1) * A(1,2) * A(2,3) * A(3,0) // 1 2 3 0 - 
    //            + A(0,2) * A(1,1) * A(2,3) * A(3,0) // 2 1 3 0 + 
    //            - A(0,2) * A(1,0) * A(2,3) * A(3,1) // 2 0 3 1 - 
               
    //            + A(0,0) * A(1,3) * A(2,1) * A(3,2) // 0 3 1 2 + 
    //            - A(0,0) * A(1,3) * A(2,2) * A(3,1) // 0 3 2 1 - 
    //            - A(0,1) * A(1,3) * A(2,0) * A(3,2) // 1 3 0 2 - 
    //            + A(0,1) * A(1,3) * A(2,2) * A(3,0) // 1 3 2 0 + 
    //            - A(0,2) * A(1,3) * A(2,1) * A(3,0) // 2 3 1 0 - 
    //            + A(0,2) * A(1,3) * A(2,0) * A(3,1) // 2 3 0 1 + 
               
    //            - A(0,3) * A(1,0) * A(2,1) * A(3,2) // 3 0 1 2 - 
    //            + A(0,3) * A(1,0) * A(2,2) * A(3,1) // 3 0 2 1 + 
    //            + A(0,3) * A(1,1) * A(2,0) * A(3,2) // 3 1 0 2 + 
    //            - A(0,3) * A(1,1) * A(2,2) * A(3,0) // 3 1 2 0 - 
    //            + A(0,3) * A(1,2) * A(2,1) * A(3,0) // 3 2 1 0 + 
    //            - A(0,3) * A(1,2) * A(2,0) * A(3,1) // 3 2 0 1 - 
    //            ;
        
    } else if( 2 <= A.getdimin() and A.getdimin() <= 8 ) {
      
        return Determinant_laplaceexpansion( A );
        
    } else {
        
        return Determinant_bareiss( A );
        
    }
}





Float Determinant_laplaceexpansion( const DenseMatrix& A )
{
    assert( A.issquare() );
    
    if( A.getdimin() == 0 )
        return 1.;
    
    Float ret = 0.;
    int signum  = 1;
    
    int i = 77;
    std::vector<int>  aux( A.getdimin() );
    std::vector<int> perm( A.getdimin() );
    for( int j = 0; j < perm.size(); j++ ) perm[j] = j;
    
    HeapsAlgorithmInit( i, aux, perm );
    
    do {
        
        Float summand = signum;
        for( int r = 0; r < A.getdimout(); r++ )
            summand *= A( r, perm[r] );
        
        ret += summand;
        
        signum = -signum;
        
    } while ( HeapsAlgorithmStep( i, aux, perm ) );
    
    return ret;
    
}



Float Determinant_gauss( DenseMatrix A )
{
    assert( A.issquare() );
    
    if( A.getdimin() == 0 )
        return 1.;
    
    const int n = A.getdimin();
    
    int signum = 1;
    
    for( int i = 0; i < n; i++ )
    {
        
        int r = i, c = i;
        for( int s = i; s < n; s++ )
        for( int d = i; d < n; d++ )
            if( absolute(A(s,d)) > absolute(A(r,c)) ) {
                r = s; c = d;
            }
        
        // make swappings in the range i..n
        assert( r >= i and r >= i );
        if( r != i ) { signum = -signum; A.swaprow   ( r, i ); }
        if( c != i ) { signum = -signum; A.swapcolumn( c, i ); }
        
        for( int s = i; s < n; s++ )
        for( int d = i; d < n; d++ )
            assert( absolute(A(s,d)) <= absolute(A(i,i)) );
        
        if( absolute( A(i,i) ) == 0.0 )
            return 0.; 

        for( int k = i+1; k < n; k++ ) {
            
            Float coeff = - A(k,i) / A(i,i);
            
            for( int j = i+1; j < n; j++ )
                A( k, j ) = A( k, j ) + coeff * A( i, j );
            
        }
    }
    
    Float ret = signum;
    for( int i = 0; i < n; i++ ) ret *= A(i,i);
    
    return ret;
    
}





Float Determinant_bareiss( DenseMatrix A )
{
    assert( A.issquare() );
    
    if( A.getdimin() == 0 )
        return 1.;
    
    const int n = A.getdimin();
    
    
    int signum = 1;

    for( int i = 0; i < n - 1; i++ )
    {

        int r = i, c = i;

        for( int s = i; s < n; s++ )
        for( int d = i; d < n; d++ )
            if( absolute(A(s,d)) > absolute(A(r,c)) ) {
                r = s; c = d;
            }
        
        // make swappings in the range i..n
        assert( r >= i and r >= i );
        if( r != i ) { signum = -signum; A.swaprow   ( r, i ); }
        if( c != i ) { signum = -signum; A.swapcolumn( c, i ); }
        
        for( int s = i; s < n; s++ )
        for( int d = i; d < n; d++ )
            assert( absolute(A(s,d)) <= absolute(A(i,i)) );
        
        if( absolute( A(i,i) ) == 0.0 )
            return 0.; 

        

        //Apply formula
        for (int r = i + 1; r < n; r++) 
        for (int c = i + 1; c < n; c++) 
        {
            A(r,c) = A(i,i) * A(r,c) - A(r,i) * A(i,c);
            if(i != 0) { A(r,c) /= A(i-1,i-1); }
        }
    }

    return signum * A(n-1,n-1);
}













/************************************
*****                          ******
*****   COFACTOR MATRIX        ******
*****                          ******
************************************/

DenseMatrix CofactorMatrix( const DenseMatrix& A )
{
  assert( A.issquare() );
  
  if( A.getdimin() == 0 ) 
    return DenseMatrix( 0, 0 );
  
  int iter = 77;
  std::vector<int>  aux( A.getdimin() - 1 );
  std::vector<int> perm( A.getdimin() - 1 );
  for( int j = 0; j < perm.size(); j++ ) perm[j] = j;
  
  DenseMatrix cof( A.getdimin(), A.getdimin(), 0. );
  
  HeapsAlgorithmInit( iter, aux, perm );
  
  int sign_perm = 1;
  
  do {
    
    for( int r = 0; r < A.getdimin(); r++ )
    for( int c = 0; c < A.getdimin(); c++ )
    {
      
      int sign_entry = sign_power( r+c );
      
      assert( sign_perm * sign_entry == 1 or sign_perm * sign_entry == -1 );
      Float summand = sign_perm * sign_entry;
      
      for( int j = 0; j < A.getdimin() - 1; j++ )
        summand *= A( j < c ? j : j+1, perm[j] < r ? perm[j] : perm[j]+1 );
      
      cof(r,c) += summand;
      
    }
    
    sign_perm *= -1;
    
  } while ( HeapsAlgorithmStep( iter, aux, perm ) );
  
  return cof;
}







/************************************
*****                          ******
*****    INVERSE               ******
*****                          ******
************************************/



DenseMatrix Inverse( DenseMatrix A )
{
    
    assert( A.issquare() );
    if( A.getdimin() == 1 ) {

        return DenseMatrix( 1, 1, 1. / A(0,0) );
    
    // } else if( A.getdimin() == 2 ) {

    //     Float D = A(0,0) * A(1,1) - A(0,1) * A(1,0);
        
    //     return DenseMatrix( 2, 2, { A(1,1)/D, -A(0,1)/D, -A(1,0)/D, A(0,0)/D } );

    // } else if( A.getdimin() == 3 ) {

    //     Float D_00 = A(1,1) * A(2,2) - A(2,1) * A(1,2);
    //     Float D_01 = A(1,0) * A(2,2) - A(2,0) * A(1,2);
    //     Float D_02 = A(1,0) * A(2,1) - A(2,0) * A(1,1);
        
    //     Float D_10 = A(0,1) * A(2,2) - A(2,1) * A(0,2);
    //     Float D_11 = A(0,0) * A(2,2) - A(2,0) * A(0,2);
    //     Float D_12 = A(0,0) * A(2,1) - A(2,0) * A(0,1);
        
    //     Float D_20 = A(0,1) * A(1,2) - A(1,1) * A(0,2);
    //     Float D_21 = A(0,0) * A(1,2) - A(1,0) * A(0,2);
    //     Float D_22 = A(0,0) * A(1,1) - A(1,0) * A(0,1);

    //     Float D = A(0,0) * D_00 - A(0,1) * D_01 + A(0,2) * D_02;
        
    //     return DenseMatrix( 3, 3, {
    //        +D_00/D, -D_10/D, +D_20/D,
    //        -D_01/D, +D_11/D, -D_21/D,
    //        +D_02/D, -D_12/D, +D_22/D           
    //     } );

    // } else if( A.getdimin() <= 6 ) {
    //     Inverse_CramersRule_InSitu( A );
    } else {
        Inverse_gauss_InSitu( A );
    }
    return A;
}


void Inverse_InSitu( DenseMatrix& A )
{
    Inverse_gauss_InSitu( A );
}





void Inverse_CramersRule_InSitu( DenseMatrix& A )
{
    assert( A.issquare() ); 
    // Float det = Determinant_laplaceexpansion( A );
    // assert( absolute( det ) > machine_epsilon );
    
    // auto B = CofactorMatrix( A ); B /= Determinant( A );
    // auto C = ( CofactorMatrix( A ) /  Determinant( A ) );
    // LOG << A << nl << B << nl << C << nl;
    // LOG << A*B << nl << A*C << nl;
    // assert( ( C - B ).is_numerically_small() );

    A = CofactorMatrix( A ) /  Determinant( A );
}


 
void Inverse_gauss_InSitu( DenseMatrix& mat, bool pivoting )
{
    
    assert( mat.issquare() );
    
    const int n = mat.getdimout();
    
    int* pivots = nullptr;
    if(pivoting) pivots = new (std::nothrow) int[n];
    
    for( int i = 0; i < n; i++ ) {
        
        if( pivoting ) {
            
            int c_max = i;
            for( int c = i+1; c < n; c++ )
                if( absolute(mat(i,c)) > absolute(mat(i,c_max)) ) 
                    c_max = c;
            
            pivots[i] = c_max;
            mat.swapcolumn( c_max, i );
            
        }
        
        for( int k = 0; k < n; k++ ) { // each
            
            if( i == k ) continue; 
            
            assert( absolute(mat(i,i)) != 0.0 );
            
            Float coeff = - mat( k, i ) / mat( i, i );
            
            for( int j = i+1; j < n; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
            }

            mat( k, i ) = coeff;
            
            for( int j = 0; j < i; j++ ) {
                mat( k, j ) = mat( k, j ) + coeff * mat( i, j );
            }
            
        }
        
        Float coeff = 1. / mat(i,i);
        
        for( int j = 0; j < i; j++ ) {
            mat(i,j) *= coeff;
        }
        
        mat(i,i) = coeff;
        
        for( int j = i+1; j < n; j++ ) {
            mat(i,j) *= coeff;
        }
            
    }
    
    if( pivoting ) {
        for( int i = n-1; i >= 0; i-- )
//         for( int i = 0; i < n; i++ ) 
        {
//             LOG << "swap " << i << space << pivots[i] << nl;
            mat.swaprow( i, pivots[i] );
        }
    }
    
    // finished!
    
    if( pivoting ) delete[] pivots;
    
}












// void InverseAndDeterminant( const DenseMatrix& A, DenseMatrix& Ainv, Float& Adet )
// {
//     assert( A.issquare() ); 
//     DenseMatrix Cof = CofactorMatrix( A );
//     Adet = Determinant( A );
//     Ainv = Cof / Adet;
//     
// }






/************************************
*****                          ******
*****  SUBDETERMINANT MATRIX   ******
*****                          ******
************************************/


DenseMatrix SubdeterminantMatrixSquare( const DenseMatrix& A, int k )
{
    A.check();
    assert( A.issquare() );
    assert( 0 <= k && k <= A.getdimin() );

    // performance "hacks", which can be disabled at any time
    if( k == 0 ) return DenseMatrix(1,1,1.);
    if( k == 1 ) return A;
    
    const int n = A.getdimin();
    IndexRange fromrange = IndexRange( 0, k-1 );
    IndexRange torange = IndexRange( 0, n-1 );
    std::vector<IndexMap> sigmas = generateSigmas( fromrange, torange );
    
    DenseMatrix ret( SIZECAST( sigmas.size() ) );
    for( int rim = 0; rim < sigmas.size(); rim++ )
    for( int cim = 0; cim < sigmas.size(); cim++ )
    {
        ret(rim,cim) = Determinant( A.submatrix( sigmas.at(rim), sigmas.at(cim) ) );
    }
    
    ret.check();
    return ret;
}


DenseMatrix SubdeterminantMatrix( const DenseMatrix& A, int k )
{
    A.check();
    assert( 0 <= k && k <= A.getdimin() && k <= A.getdimout() );

    if( A.issquare() ) return SubdeterminantMatrixSquare(A,k);

    // performance "hacks", which can be disabled at any time
    if( k == 0 ) return DenseMatrix(1,1,1.);
    if( k == 1 ) return A;

    const IndexRange range_from = IndexRange( 0, k-1 );
    const IndexRange range_rows = IndexRange( 0, A.getdimout()-1 );
    const IndexRange range_cols = IndexRange( 0, A.getdimin()-1 );
    const std::vector<IndexMap> sigmas_rows = generateSigmas( range_from, range_rows );
    const std::vector<IndexMap> sigmas_cols = generateSigmas( range_from, range_cols );
    
    DenseMatrix ret( SIZECAST( sigmas_rows.size() ), SIZECAST( sigmas_cols.size() ), 0. );
    for( int rim = 0; rim < sigmas_rows.size(); rim++ )
    for( int cim = 0; cim < sigmas_cols.size(); cim++ )
    {
        ret(rim,cim) = Determinant( A.submatrix( sigmas_rows.at(rim), sigmas_cols.at(cim) ) );
    }
    
    ret.check();
    return ret;
}






DenseMatrix MatrixTensorProduct( const DenseMatrix& A, const DenseMatrix& B )
{
    A.check(); B.check();
    const int newrows = A.getdimout() * B.getdimout();
    const int newcols = A.getdimin()  * B.getdimin();
    
    DenseMatrix ret( newrows, newcols );
    assert( ret.getdimout() == newrows );
    assert( ret.getdimin()  == newcols );

    const int num_A_rows = A.getdimout();
    const int num_A_cols  = A.getdimin();
    const int num_B_rows = B.getdimout();
    const int num_B_cols  = B.getdimin();
    
    for( int rl = 0; rl < num_A_rows; rl++ )
    for( int cl = 0; cl < num_A_cols; cl++ )
    for( int rr = 0; rr < num_B_rows; rr++ )
    for( int cr = 0; cr < num_B_cols; cr++ )
    {
        ret( rl * num_B_rows + rr, cl * num_B_cols + cr )
        =
        A( rl, cl ) * B( rr, cr );
    }
    
    ret.check();
    return ret;
}

