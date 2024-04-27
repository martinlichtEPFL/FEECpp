#ifndef INCLUDEGUARD_FEM_TRACEMATRIX
#define INCLUDEGUARD_FEM_TRACEMATRIX


#include "../sparse/sparsematrix.hpp"
#include "../mesh/mesh.hpp"


SparseMatrix FEECBrokenTraceMatrix( const Mesh& mesh, int n, int k, int r, bool is_signed );



#endif
