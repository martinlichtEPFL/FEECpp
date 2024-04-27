#ifndef INCLUDEGUARD_FEM_FEECSULLIVANINCLUSIONMATRIX
#define INCLUDEGUARD_FEM_FEECSULLIVANINCLUSIONMATRIX

#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous Sullivan forms                       //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECSullivanInclusionMatrix( const Mesh& mesh, int n, int k, int r );




#endif
