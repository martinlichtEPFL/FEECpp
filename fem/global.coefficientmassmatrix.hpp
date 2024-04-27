#ifndef INCLUDEGUARD_FEM_FEECBROKENCOEFFICIENTMASSMATRIX
#define INCLUDEGUARD_FEM_FEECBROKENCOEFFICIENTMASSMATRIX


#include <functional>

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

SparseMatrix FEECBrokenCoefficientMassMatrix( const Mesh& mesh, int n, int k, int r, 
                                              int w, const std::function<DenseMatrix(const FloatVector&)>& generator );

#endif
