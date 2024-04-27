#ifndef INCLUDEGUARD_FEM_LAGRANGE_MATRICES
#define INCLUDEGUARD_FEM_LAGRANGE_MATRICES

#include <functional>

#include "../sparse/sparsematrix.hpp"
#include "../sparse/matcsr.hpp"
#include "../mesh/mesh.hpp"






//////////////////////////////////////////////////////
//                                                  //
//  Matrix for the broken Lagrange mass pairing     //
//  with poly degree r                              //
//                                                  //
//////////////////////////////////////////////////////

SparseMatrix LagrangeBrokenMassMatrix( const Mesh& mesh, int r );



///////////////////////////////////////////////////////
//                                                   //
//  Matrix for the continuous Lagrange mass pairing  //
//  with poly degree r                               //
//                                                   //
///////////////////////////////////////////////////////

SparseMatrix LagrangeMassMatrix( const Mesh& mesh, int r );
MatrixCSR LagrangeCoefficientMassMatrix( const Mesh& mesh, int r, int w, const std::function<Float(const FloatVector&)> weight );



//////////////////////////////////////////////////////////
//                                                      //
//  Matrix for the broken Lagrange stiffness pairing    //
//  with poly degree r                                  //
//                                                      //
//////////////////////////////////////////////////////////

SparseMatrix LagrangeBrokenStiffnessMatrix( const Mesh& mesh, int r );






//////////////////////////////////////////////////////////////
//                                                          //
//  Matrix for the continuous Lagrange stiffness pairing    //
//  with poly degree r                                      //
//                                                          //
//////////////////////////////////////////////////////////////

SparseMatrix LagrangeStiffnessMatrix( const Mesh& mesh, int r );
MatrixCSR LagrangeCoefficientStiffnessMatrix( const Mesh& mesh, int r, int w, const std::function<DenseMatrix(const FloatVector&)> weight );






////////////////////////////////////////////////
//                                            //
//  Matrix for the inclusion of               //
//  continuous Lagrange elements              //
//  into the broken space                     //
//  with poly degree r                        //
//                                            //
////////////////////////////////////////////////

SparseMatrix LagrangeInclusionMatrix( const Mesh& mesh, int n, int r );







#endif
