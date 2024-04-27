#ifndef INCLUDEGUARD_FEM_FEECBROKENELEVATIONMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENELEVATIONMATRIX


#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for degree elevantion                    //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//  and lifts up rplus degrees                      //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECBrokenElevationMatrix( const Mesh& mesh, int n, int k, int r, int r_plus );

#endif
