#ifndef INCLUDEGUARD_FEM_FEECCECHMATRIX
#define INCLUDEGUARD_FEM_FEECCECHMATRIX


#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"




//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the mass pairing                     //
//  on broken spaces of P_r \Lambda^k               //
//                                                  //
//  Works on any n dimensional mesh                 //
//  with any form degree k and poly degree r        //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix FEECCechMassMatrix( const Mesh& mesh, int n, int k, int s );

SparseMatrix FEECCechDiffMatrix( const Mesh& mesh, int n, int k );


#endif
