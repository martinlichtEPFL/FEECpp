
#include <algorithm> // for shuffle 
#include <random> // for shuffle 
#include <vector>

#include "../basic.hpp"
#include "matcsr.hpp"

#include "rainbow.hpp"



Rainbow::Rainbow( const MatrixCSR& mat, bool do_shuffle )
{
    mat.check();
    assert( mat.issquare() );

    // const std::vector<int>& A = mat.getA();
    // const std::vector<int>& C = mat.getC();
    const int* __restrict__ A = mat.getA();
    const int* __restrict__ C = mat.getC();

    const int N = mat.getdimout();

    num_rows = N;

    
    // 1. Determine the number of colors and assign each row its color (F)

    num_colors = 1;

    std::vector<int> painting_order( N ); // determine in which orders the colors are assigned

    for( int i = 0; i < N; i++ ) painting_order[i] = i;
    
    if( do_shuffle) {
        std::random_device rd;
        std::mt19937 g(rd());
        std::shuffle( painting_order.begin(), painting_order.end(), g );
    }

    F = std::vector<int>( N, 0 );

    for( int p = 0; p < N; p++ )
    {
        
        int r = painting_order[p];
        
        loop:
        for( int c = A[r]; c < A[r+1]; c++ )
        {
            int j = C[c];
            
            if( j == r ) continue;

            if( F[j] == F[r] )
            {
                F[r]++;
                goto loop;
            }
        }

        if( F[r] == num_colors ) num_colors++;
        
        assert( 0 <= F[r] and F[r] < num_colors );
        for( int c = A[r]; c < A[r+1]; c++ )
        {
            int j = C[c];
            if( j != r ) assert( F[j] != F[r] );
        }

    }

    for( int r = 0; r < N; r++ ) assert( 0 <= F[r] and F[r] < num_colors );



    // 2. count how many rows there are of each color 

    std::vector<int> counter( num_colors );

    for( int r = 0; r < N; r++ ) counter[ F[r] ]++;
    
    for( int c = 0; c < num_colors; c++ ) assert( 0 < counter[c] and counter[c] <= N );

    

    // 3. create the array of limits (B)

    B = std::vector<int>( num_colors+1, 0 );

    for( int c = 1; c < num_colors+1; c++ ) B[c]  = counter[c-1];
    for( int c = 1; c < num_colors+1; c++ ) B[c] += B[c-1];

    Assert( B[num_colors] == N, num_colors, B[num_colors], N );
    for( int c = 0; c < num_colors; c++ ) assert( B[c+1] == counter[c] + B[c] );
    

    // 4. fill in the row indices (R)

    R = std::vector<int>( N, -1 );

    for( int r = 0; r < N; r++ )
    {
        int color = F[r];
        assert( 0 <= color and color < num_colors );
        assert( counter[color] > 0 );
        counter[color]--;
        assert( B[color] + counter[color] < B[color+1] );
        assert( B[color] + counter[color] < N );
        assert( R[ B[color] + counter[color] ] == -1 );
        R[ B[color] + counter[color] ] = r;
    }

    for( int c = 0; c < num_colors; c++ ) assert( counter[c] == 0 );


    // Assert that the colors are all the correct ones 

    check();

    // for( int t : B ) LOG << t << nl;
    // LOG << "number rows: " << N << nl;
    // LOG << "number colors: " << num_colors << nl;

}


void Rainbow::check() const
{
    const int N = num_rows;

    for( int c = 0; c < num_colors; c++ )
        for( int i = B[c]; i < B[c+1]; i++ )
            assert( F[ R[ i ] ] == c );

    for( int r = 0; r < N; r++ )
    {
        bool found = false;
        
        int color = F[r];

        for( int s = B[color]; not found and s < B[color+1]; s++ )
            if( R[s] == r ) found = true;

        assert( found );
    }
}

std::string Rainbow::text() const 
{
    std::string str;
    str += "num_colors: ";
    str += std::to_string(num_colors);
    return str;
}
