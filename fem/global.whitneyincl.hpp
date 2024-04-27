#ifndef INCLUDEGUARD_FEM_FEECWHITNEYINCLUSIONMATRIX
#define INCLUDEGUARD_FEM_FEECWHITNEYINCLUSIONMATRIX

#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous Whitney forms                        //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //  
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECWhitneyInclusionMatrix( const Mesh& mesh, int n, int k, int r );




#endif
