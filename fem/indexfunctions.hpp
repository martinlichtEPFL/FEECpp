#ifndef INCLUDEGUARD_FEM_INDEXFUNCTIONS_HPP
#define INCLUDEGUARD_FEM_INDEXFUNCTIONS_HPP


#include <utility>
#include <vector>

#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"



//////////////////////////////////////////////////////
//                                                  //
//  List of Sullivan Indices                        //
//                                                  //
//  Constructs the basis for the Sullivan space     //
//                                                  //
//  alpha is a multiindex                           //
//  sigma 1:k -> 0:n                                //
//  min[alpha] notin [sigma]                        //
//  [alpha] u [sigma] = [0..n]                      //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

std::vector< std::pair<MultiIndex,IndexMap> > ListOfSullivanIndices( int n, int k, int r );

//////////////////////////////////////////////////////
//                                                  //
//  List of Sullivan Indices                        //
//                                                  //
//  Constructs the basis for the Whitney space      //
//                                                  //
//  alpha is a multiindex                           //
//  rho   0:k -> 0:n                                //
//  0 in [rho]                                      //
//  [alpha] u [rho] = [0..n]                        //
//                                                  //
//                                                  //
//////////////////////////////////////////////////////

std::vector< std::pair<MultiIndex,IndexMap> > ListOfWhitneyIndices( int n, int k, int r );



#endif
