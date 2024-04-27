#ifndef INCLUDEGUARD_SORTHACK_HPP
#define INCLUDEGUARD_SORTHACK_HPP

#include <cassert>

#include <algorithm>
#include <utility>
#include <vector>

#include <vector>

template<typename T>
void insertionSort( std::vector<T>& vec, int low, int high )
{
    assert( 0 <= low and low <= high and high <= vec.size() - 1 );

    for( int i = low + 1; i <= high; i++ )
    {
        T key = vec[i];
        int j = i - 1;

        while( j >= low && vec[j] > key )
        {
            vec[j + 1] = vec[j];
            j = j - 1;
        }
        vec[j + 1] = key;
    }

    for( int i = low + 1; i <= high; i++ )
        assert( vec[i - 1] <= vec[i] );
}

/*
template<typename T>
void insertionSort( std::vector<T>& vec, int low, int high )
{
    assert( 0 <= low and low <= high and high <= vec.size() - 1 );

    for( int i = low + 1; i <= high; i++ )
    {
        T key = vec[i];
        
        int j = i;

        while( j > low && vec[j-1] > key )
        {
            vec[j] = vec[j-1];
            j--;
        }
        assert( j >= low );
        assert( vec[j] >= key );
        vec[j] = key;
    }

    for( int i = low + 1; i <= high; i++ )
        assert( vec[i - 1] <= vec[i] );
}
*/

template<typename T>
void myswap( T* a, T* b )
{
    T t = *a;
    *a = *b;
    *b = t;
}

template<typename T>
int partitionCGPT( std::vector<T>& vec, int low, int high )
{
    assert( 0 <= low and low <= high and high <= vec.size() - 1 );

    T pivot = vec[low];
    int i = low - 1, j = high + 1;

    while( true )
    {
        do {
            i++;
        } while( vec[i] < pivot );

        do {
            j--;
        } while( vec[j] > pivot );

        if( i >= j )
            return j;

        myswap( &vec[i], &vec[j] );
    }
}




template<typename T>
void quickSortCGPT( std::vector<T>& vec, int low, int high )
{
    assert( 0 <= low and low <= high and high <= vec.size() - 1 );

    if( low < high )
    {

        if( high - low <= 32 )
        {
            insertionSort( vec, low, high );
            return;
        }

        int pi = partitionCGPT( vec, low, high );

        // assert( low <= pi && pi+1 <= high );
        for( int i =  low; i <=   pi; i++ ) 
        for( int j = pi+1; j <= high; j++ ) 
           assert( vec[i] <= vec[j] );


        quickSortCGPT( vec,    low,   pi );
        quickSortCGPT( vec, pi + 1, high );

        for( int i = low + 1; i <= high; i++ )
            assert( vec[i - 1] <= vec[i] );
    }
}

// QuickSort function
template<typename T>
void quickSortCGPT( std::vector<T>& vec )
{
    quickSortCGPT( vec, 0, vec.size() - 1 );
}

// template<typename T>
// void printArray( std::vector<T>& vec, int size )
// {
//     for( int i = 0; i < size; i++ ) std::cout << vec[i] << " ";
//     std::cout << nl;
// }

// // Driver code
// int main() {
//     std::vector<int> vec = { 10, 7, 8, 9, 1, 5 };
//     int n = vec.size();
//     quickSort( vec, 0, n - 1 );
//     std::cout << "Sorted vecay: \n";
//     printArray( vec, n );
//     return 0;
// }

template<typename A>
inline void sorthack( std::vector<A>& vec )
{

    // quickSort( vec );
    std::sort( vec.begin(), vec.end() );
    /*
    //     for( int s = 0; s < vec.size(); s++ )
    //     for( int t = 0; t < vec.size(); t++ )
    //     {
    //         if( s < t && vec[s] > vec[t] )
    //         {
    //             A temp = vec[t];
    //             vec[t] = vec[s];
    //             vec[s] = temp;
    //         }
    //     }

        // insertion sort
        for( int s = 0; s < vec.size(); s++ )
        for( int t = s; t > 0 && vec.at(t-1) > vec.at(t); t-- )
        {
            A temp = vec[t-1];
            vec[t-1] = vec[t];
            vec[t] = temp;
        }

        for( int t = 1; t < vec.size(); t++ )
            assert( vec[t-1] <= vec[t] );
    */
}

#endif
