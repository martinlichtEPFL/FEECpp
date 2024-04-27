#ifndef INCLUDEGUARD_SOLVER_NULLSPACE
#define INCLUDEGUARD_SOLVER_NULLSPACE

#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"


std::vector<FloatVector> computeNullspace(
    const LinearOperator& SystemMatrix,
    const LinearOperator& MassMatrix,
    const int max_number_of_candidates,
    //
    const Float tolerance_residual, 
    const Float tolerance_zero
);

#endif
