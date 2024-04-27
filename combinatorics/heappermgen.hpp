#ifndef INCLUDEGUARD_COMBINATORICS_HEAPPERMGEN_HPP
#define INCLUDEGUARD_COMBINATORICS_HEAPPERMGEN_HPP


#include <vector>

#include "../basic.hpp"


/***************
*** 
***  Generate Permutations using Heaps Algorithm
***  
***  The variables `seed` and `memo` are just auxiliary data,
***  whereas `perm` contains the proper data.
***  
***  Init - initiliazes the auxiliary data: `seed` and `memo`
***  
***  Step - Returns true if a has experienced a transposition
***         Returns false if a has not been changed.
***         Supposed to be used in a do-while loop.
***  
***  
***  
***  Typicall usage:
***  
***  HeapsAlgorithmInit( seed, memo, perm );
***  
***  do { 
***  
***     do_stuff_with( perm );
***  
***  } while( HeapsAlgorithmInit( seed, memo, perm ) );
***  
***  
***************/

void HeapsAlgorithmInit( int& seed, std::vector<int>& memo, const std::vector<int>& perm );

bool HeapsAlgorithmStep( int& seed, std::vector<int>& memo,       std::vector<int>& perm );


struct HeapAlgorithmState
{
    int seed;
    std::vector<int> memo;
    std::vector<int> perm;
};

inline void HeapsAlgorithmInit( HeapAlgorithmState& state )
{
    HeapsAlgorithmInit( state.seed, state.memo, state.perm );
}

inline bool HeapsAlgorithmStep( HeapAlgorithmState& state )
{
    return HeapsAlgorithmStep( state.seed, state.memo, state.perm );
}





#endif
