#ifndef INCLUDEGUARD_COMBINATORICS_GENERATEINDEXMAPS_HPP
#define INCLUDEGUARD_COMBINATORICS_GENERATEINDEXMAPS_HPP


#include <vector>

#include "../basic.hpp"

#include "indexrange.hpp"
#include "indexmap.hpp"

/***************
*** 
***  Generate Index Maps of different kinds 
*** 
***  0) generate empty mapping 
***  1) generate all possible mappings 
***  2) generate Permuations, and tell their sign 
***  3) generate the Sigma mappings 
***  
***************/





std::vector<IndexMap> generateEmptyMap( const IndexRange& from, const IndexRange& to );


/*****
 * 
 * Generate all the indexmaps from range `from` mapping into the range `to`.
 * 
 * If `from` is empty, then there exists exactly one such map,
 * namely the empty map.
 * 
 * If `from` is not empty but `to` is empty, 
 * then no such map exists (for no input value can we fix an output value).
 * 
 * If none of the two ranges are empty, the output is the obvious one. 
 * The cardinality of the index mappings is 
 *      [ size of output range ] ** [ size of input range ]
 * 
 * NOTE:
 * In the special case where both ranges are empty,
 * we'd set 0 ** 0 := 0.
 * 
 ****/

std::vector<IndexMap> generateIndexMaps( const IndexRange& range );
std::vector<IndexMap> generateIndexMaps( const IndexRange& from, const IndexRange& to );


/*****
 * 
 * Generate all permutations of the range `ir`.
 * 
 * If `ir` is empty, then there exists exactly one such permutation,
 * namely the empty one. Its sign is 1 just like the determinant 
 * of the empty matrix is 1.
 * 
 ****/

std::vector<IndexMap> generatePermutations( const IndexRange& ir );

int signPermutation( const IndexMap& im );


/*****
 * 
 * Generate all the -ascending- indexmaps from range `from` mapping into the range `to`.
 * 
 * If `from` is empty, then there exists exactly one such map, namely the empty map.
 * 
 * If `from` is not empty but `to` is empty, then no such map exists 
 * (for no input value can we fix an output value).
 * 
 * If none of the two ranges are empty, the output is the obvious one. 
 * The cardinality of the index mappings is 
 *      [ size of output range ] \choose [ size of input range ]
 * 
 * NOTE:
 * We always have ( n \choose 0 ) = 1 regardless of (non-negative) n.
 * We always have ( n \choose k ) = 0 if k > n and n non-negative.
 * 
 * 
 * 
 ****/

std::vector<IndexMap> generateSigmas( const IndexRange& from, const IndexRange& to );



#endif
