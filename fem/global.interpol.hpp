#ifndef INCLUDEGUARD_FEM_FEECBROKENINTERPOLATIONMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENINTERPOLATIONMATRIX


#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for degree reduction                     //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//  and reduces by rplus degrees                    //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECBrokenInterpolationMatrix( const Mesh& mesh, int n, int k, int r, int r_plus );

#endif
