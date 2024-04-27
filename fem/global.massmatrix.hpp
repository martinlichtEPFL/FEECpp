#ifndef INCLUDEGUARD_FEM_FEECBROKENMASSMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENMASSMATRIX


#include "../operators/floatvector.hpp"
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

SparseMatrix FEECBrokenMassMatrix( const Mesh& mesh, int n, int k, int r );

SparseMatrix FEECBrokenMassMatrixRightFactor( const Mesh& mesh, int n, int k, int r );

FloatVector FEECBrokenMassMatrix_cellwisemass( const Mesh& mesh, int n, int k, int r, const FloatVector vec );

SparseMatrix FEECBrokenMassMatrix_cellwiseinverse( const Mesh& mesh, int n, int k, int r );


#endif
