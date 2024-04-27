#ifndef INCLUDEGUARD_FEM_FEECBROKENDIFFMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENDIFFMATRIX

#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the exterior derivative              //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECBrokenDiffMatrix( const Mesh& mesh, int n, int k, int r );




#endif
