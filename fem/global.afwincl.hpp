#ifndef INCLUDEGUARD_FEM_FEECAFWINCLUSIONMATRIX
#define INCLUDEGUARD_FEM_FEECAFWINCLUSIONMATRIX

#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


//////////////////////////////////////////////////////
//                                                  //
//  Matrix for inclusion of                         //
//  continuous AFW forms                            //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //  
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECAFWInclusionMatrix( const Mesh& mesh, int n, int k, int r );




#endif
