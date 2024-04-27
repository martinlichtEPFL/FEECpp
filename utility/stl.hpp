#ifndef INCLUDEGUARD_UTILITY_STL_HPP
#define INCLUDEGUARD_UTILITY_STL_HPP

#include <algorithm>
#include <functional>
#include <vector>

#include "../basic.hpp"


/////////////////////////////////////////////////
//                                             //
//       SUM INTEGERS PRODUCED BY LAMBDA       //
//                                             //
/////////////////////////////////////////////////

inline int sum_int( int from, int to, const std::function<int(int)>& calc )
{
    if( from > to )
        return 0;
    int ret = 0;
    for( int i = from; i <= to; i++ )
        ret += calc( i );
    return ret;
}

inline int sum_int( int to, const std::function<int(int)>& calc )
{
    return sum_int( 0, to, calc );
}




/******************************************************/
/*        write zero-based range into vector          */
/******************************************************/

inline std::vector<int> range( int to )
{
    Assert( to >= 0 );
    std::vector<int> ret(to+1);
    for( int i = 0; i <= to; i++ ) ret.at(i) = i;
    Assert( ret.size() == to+1 );
    return ret;
}


/******************************************************/
/*   remove duplicates from random access container   */
/******************************************************/

template< typename T >
inline void sort_and_remove_duplicates( T& t )
{
    std::sort( t.begin(), t.end() );
    auto last = std::unique( t.begin(), t.end() );
    t.erase( last, t.end() );
}


/******************************************************/
/*      find index of element with STL vector         */
/******************************************************/

template<typename T>
inline int find_index( const std::vector<T>& vec, const T& t )
{
   const auto it = std::find( vec.begin(), vec.end(), t );
   Assert( it != vec.end() );
   const auto ret = std::distance( vec.begin(), it );
   Assert( ret >= 0 );
   Assert( ret < vec.size() );
   return SIZECAST( ret );
}


/******************************************************/
/*            merge two sorted STL lists              */
/******************************************************/

// template<typename T>
// inline void mergeelementsinsortedlist
// ( std::list<T>& L, 
//   const std::function<T( const T&, const T& )>& merge,
//   const std::function<bool( const T&, const T& )>& compare
// ) {
//     typename std::list<T>::iterator it = L.begin();
//     while( it != L.end() ){

//         typename std::list<T>::iterator now = it; 
//         typename std::list<T>::iterator next = ++it;

//         if( next == L.end() ) return;

//         if( is_equal_to( *it, *next ) )
//         {
//             *now = merge( *now, *next );
//             L.erase( next );
//             it = now;
//         } 

//     }
// }

/***********************************************/
/*         make_unique HACK                    */ 
/***********************************************/

#if __cplusplus < 201402L

/****
 * 
 * A very imperfect solution for make_unique in C++11
 * We enter undefined behavior territory here
 * 
 ****/
#warning \
This code extends the std namespace so that `make_unique` is available throughout the code. \
This was triggered by a C++ version below C++14. \
The lack of `make_unique` is generally considered an oversight in the official standard. \
While this may be a practical workaround, it is officially considered undefined behavior. \
Please try to compile with C++14 or higher.

#include <memory>

namespace std
{
template <typename T, typename ...Args> 
inline std::unique_ptr<T> make_unique(Args && ...args)
{
  return std::unique_ptr<T>( new T(std::forward<Args>(args)...) );
}
}

#endif














#endif
