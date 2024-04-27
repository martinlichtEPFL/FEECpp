#ifndef INCLUDEGUARD_FEM_FEECFLAGMATRIX
#define INCLUDEGUARD_FEM_FEECFLAGMATRIX

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

DiagonalOperator FEECSullivanFlagMatrix( const Mesh& mesh, int n, int k, int r );

DiagonalOperator FEECWhitneyFlagMatrix( const Mesh& mesh, int n, int k, int r );




#endif
