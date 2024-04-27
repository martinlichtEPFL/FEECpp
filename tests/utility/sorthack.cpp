// g++ sorthacktest.cpp -o sorthacktest

#include <cassert>
#include <cstdio>
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <iostream>
#include <vector>

#include "../../utility/random.hpp"
#include "../../utility/sorthack.hpp"
#include "../../utility/profilerutils.hpp"

static void insertionsort( std::vector<int>& vec, std::size_t begin, std::size_t end )
{
    if( begin >= end )
        return;

    for( std::size_t i = begin + 1; i < end; i++ )
    {
        int thing = vec[i];
        std::size_t j = i;
        while( j > begin && vec[j - 1] > thing )
        {
            vec[j] = vec[j - 1];
            j--;
        }
        vec[j] = thing;

        for( std::size_t j = begin + 1; j <= i; j++ )
            assert( vec[j - 1] <= vec[j] );
    }

    for( std::size_t j = begin + 1; j < end; j++ )
        assert( vec[j - 1] <= vec[j] );
}

static void insertionsort( std::vector<int>& vec ) { insertionsort( vec, 0, vec.size() ); }

template<typename T>
static int partition( std::vector<T>& vec, int low, int high )
{
    assert( 0 <= low and low <= high and high <= vec.size() - 1 );

    size_t i = low;
    size_t j = high;

    T pivot = vec[high];

    for( ; 1; i++, j-- )
    {
        while( vec[i] < pivot )
            i++;
        while( vec[j] > pivot )
            j++;
        if( i >= j )
            return i;
        myswap( &vec[i], &vec[j] );
    }

    assert( false );
    while( i < j )
    {
        while( i < j && vec[i] <= pivot )
        {
            i++;
        };

        while( i < j && vec[j] > pivot )
        {
            j--;
        }

        // assert( i <= j );

        if( vec[i] > vec[j] )
            myswap( &vec[i], &vec[j] );
    }

    if( vec[i] > pivot )
        myswap( &vec[i], &vec[j] );
    else
        i = high;

    for( size_t k = low; k < i; k++ )
        assert( vec[k] <= pivot );
    for( size_t k = j; k <= high; k++ )
        assert( vec[k] >= pivot );

    return i;
}

template<typename T>
static void quickSort( std::vector<T>& vec, std::size_t low, std::size_t high )
{
    assert( 0 <= low and low <= high and high <= vec.size() - 1 );

    if( low < high )
    {

        if( high - low <= 32 )
        {
            insertionSort( vec, low, high );
            return;
        }

        std::size_t pi = partition( vec, low, high );

        // assert( low <= pi && pi+1 <= high );
        // for( int i =  low; i <=   pi; i++ )
        // for( int j = pi+1; j <= high; j++ )
        //    assert( vec[i] <= vec[j] );

        quickSort( vec, low, pi );
        quickSort( vec, pi + 1, high );

        // for( int i = low + 1; i <= high; i++ )
        //     assert( vec[i - 1] <= vec[i] );
    }
}

// QuickSort function
template<typename T>
static void quickSort( std::vector<T>& vec )
{
    quickSort( vec, 0, vec.size() - 1 );
}

// New quickSort function
template<typename T>
static void NewquickSort( std::vector<T>& vec, std::size_t low, std::size_t high )
{
    assert( low <= high );

    // printf("Length: %lu\n", (unsigned long)(high - low) );

    if( high - low < 32 ) {
        insertionsort( vec, low, high+1 );
        return;
    }

    std::size_t pivot = vec[low + ( high - low ) / 2];
    std::size_t i = low;
    std::size_t j = high;

    while( i <= j )
    {

        while( vec[i] < pivot )
            i++;
        while( vec[j] > pivot )
            j--;
        if( i <= j )
        {
            auto temp = vec[i];
            vec[i] = vec[j];
            vec[j] = temp;
            i++;
            j--;
        }

        if( low < j )
        {
            NewquickSort( vec, low, j );
        }
        if( i < high )
        {
            NewquickSort( vec, i, high );
        }
    }
}

template<typename T>
static void NewquickSort( std::vector<T>& vec )
{
    NewquickSort( vec, 0, vec.size() - 1 );
}

// Merges two subvecays of vec[].
// First subvecay is vec[l..m]
// Second subvecay is vec[m+1..r]
// Inplace Implementation
template<typename T>
void merge( std::vector<T>& vec, std::size_t start, std::size_t mid, std::size_t end )
{
    std::size_t start2 = mid + 1;

    // If the direct merge is already sorted
    if( vec[mid] <= vec[start2] )
    {
        return;
    }

    // Two pointers to maintain start
    // of both vecays to merge
    while( start <= mid && start2 <= end )
    {

        // If element 1 is in right place
        if( vec[start] <= vec[start2] )
        {
            start++;
        } else
        {
            int value = vec[start2];
            int index = start2;

            // Shift all the elements between element 1
            // element 2, right by 1.
            while( index != start )
            {
                vec[index] = vec[index - 1];
                index--;
            }
            vec[start] = value;

            // Update all the pointers
            start++;
            mid++;
            start2++;
        }
    }
}

/* l is for left index and r is right index of the
   sub-vecay of vec to be sorted */
template<typename T>
void mergeSort( std::vector<T>& vec, std::size_t l, std::size_t r )
{
    if( l < r )
    {

        // Same as (l + r) / 2, but avoids overflow
        // for large l and r
        int m = l + ( r - l ) / 2;

        // Sort first and second halves
        mergeSort( vec, l, m );
        mergeSort( vec, m + 1, r );

        merge( vec, l, m, r );
    }
}

int main( int argc, char *argv[] )
{
    std::size_t N = 2 << 10;
    std::vector<int> foo( N );

    constexpr int number_of_methods = 5;

    std::function<void( std::vector<int>& )> sorting_calls[number_of_methods] = {
        []( std::vector<int>& vec ) -> void { std::sort( vec.begin(), vec.end() ); },
        []( std::vector<int>& vec ) -> void { insertionsort( vec ); },
        []( std::vector<int>& vec ) -> void { quickSortCGPT( vec ); },
        []( std::vector<int>& vec ) -> void { mergeSort( vec, 0, vec.size() - 1 ); },
        []( std::vector<int>& vec ) -> void { NewquickSort( vec ); } };

    const char* sorting_call_names[number_of_methods] = {
        "std::sort",
        "custom insertion sort",
        "ChatGPT quick sort",
        "Merge sort",
        "Custom quick sort",
    };

    for( int c = 0; c < number_of_methods; c++ )
    {
        // skip insertion sort on long arrays
        if( c == 1 && N > 20000 ) continue;

        // skip custom quick sort on long arrays
        if( c == 4 ) continue;


        auto& sorting_call = sorting_calls[c];
        auto& name = sorting_call_names[c];

        {
            std::printf( "Stricly ascending sequence:                    %zu\t%s\n", N, name );

            for( std::size_t i = 0; i < N; i++ )
                foo[i] = i;

            {
                StopWatch watch;
                sorting_call( foo );
            }

            //         for( std::size_t i = 0; i < N; i++ ) std::cout << foo[i] << ' '; std::cout << '\n';

            for( std::size_t i = 1; i < N; i++ )
                assert( foo[i - 1] <= foo[i] );

            for( std::size_t i = 1; i < N; i++ )
                assert( foo[i] == i );
        }

        {
            std::printf( "Stricly descending sequence:                   %zu\t%s\n", N, name );

            for( std::size_t i = 0; i < N; i++ )
                foo[i] = N - i - 1;

            {
                StopWatch watch;
                sorting_call( foo );
            }

            //         for( std::size_t i = 0; i < N; i++ ) std::cout << foo[i] << ' '; std::cout << '\n';

            for( std::size_t i = 1; i < N; i++ )
                assert( foo[i - 1] <= foo[i] );

            for( std::size_t i = 1; i < N; i++ )
                assert( foo[i] == i );
        }

        {
            std::printf( "Random reorder of strictly ascending sequence: %zu\t%s\n", N, name );

            for( std::size_t i = 0; i < N; i++ )
                foo[i] = i;

            //         for( std::size_t i = 0; i < N; i++ ) std::cout << foo[i] << ' '; std::cout << '\n';

            for( std::size_t i = 0; i < N; i++ )
                std::swap( foo[i], foo[i + random_integer() % ( N - i )] );

            //         for( std::size_t i = 0; i < N; i++ ) std::cout << foo[i] << ' '; std::cout << '\n';

            {
                StopWatch watch;
                sorting_call( foo );
            }

            //         for( std::size_t i = 0; i < N; i++ ) std::cout << foo[i] << ' '; std::cout << '\n';

            for( std::size_t i = 1; i < N; i++ )
                assert( foo[i - 1] <= foo[i] );

            for( std::size_t i = 1; i < N; i++ )
                assert( foo[i] == i );
        }

        {
            std::size_t K = 1000;
            std::printf( "Randomized sequence, all values modulo %zu:   %zu\t%s\n", K, N, name );

            for( std::size_t i = 0; i < N; i++ )
                foo[i] = random_integer() % K;

            {
                StopWatch watch;
                sorting_call( foo );
            }

            for( std::size_t i = 1; i < N; i++ )
                assert( foo[i - 1] <= foo[i] ); // std::cout << foo[i] << '\n';
        }
    }

    return 0;
}
