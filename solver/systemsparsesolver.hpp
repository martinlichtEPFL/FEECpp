#ifndef INCLUDEGUARD_SOLVER_SYSTEMSPARSESOLVER
#define INCLUDEGUARD_SOLVER_SYSTEMSPARSESOLVER


#include "../basic.hpp"
#include "sparsesolver.hpp"

const Float expected_sign_of_A =  1.;

int HodgeConjugateResidualSolverCSR( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    const int*  Crows, const int*  Ccolumns, const Float*  Cvalues, 
    Float* residual,
    Float tolerance,
    int print_modulo,
    Float inneriteration_tolerance,
    int inneriteration_print_modulo
);

int HodgeConjugateResidualSolverCSR_diagonal( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    const int*  Crows, const int*  Ccolumns, const Float*  Cvalues, 
    Float* residual,
    Float tolerance,
    int print_modulo,
    Float inneriteration_tolerance,
    int inneriteration_print_modulo
);

int HodgeConjugateResidualSolverCSR_textbook( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    const int*  Crows, const int*  Ccolumns, const Float*  Cvalues, 
    Float* residual,
    Float tolerance,
    int print_modulo,
    Float inneriteration_tolerance,
    int inneriteration_print_modulo
);

int HodgeConjugateResidualSolverCSR_SSOR( 
    const int N, 
    const int L, 
    Float* x, 
    const Float* b, 
    const int*  Arows, const int*  Acolumns, const Float*  Avalues, 
    const int*  Brows, const int*  Bcolumns, const Float*  Bvalues, 
    const int* Btrows, const int* Btcolumns, const Float* Btvalues, 
    const int*  Crows, const int*  Ccolumns, const Float*  Cvalues, 
    Float* residual,
    Float tolerance,
    int print_modulo,
    Float inneriteration_tolerance,
    int inneriteration_print_modulo
);

int HodgeHerzogSoodhalterMethod( 
    const int dimension_A, 
    const int dimension_C, 
    Float* __restrict__ x_A, 
    Float* __restrict__ x_C, 
    const Float* __restrict__ b_A, 
    const Float* __restrict__ b_C, 
    const int* __restrict__  Arows, const int* __restrict__  Acolumns, const Float* __restrict__  Avalues, 
    const int* __restrict__  Brows, const int* __restrict__  Bcolumns, const Float* __restrict__  Bvalues, 
    const int* __restrict__ Btrows, const int* __restrict__ Btcolumns, const Float* __restrict__ Btvalues, 
    const int* __restrict__  Crows, const int* __restrict__  Ccolumns, const Float* __restrict__  Cvalues, 
    Float tolerance,
    int print_modulo,
    const int* __restrict__ PArows, const int* __restrict__ PAcolumns, const Float* __restrict__ PAvalues, 
    const int* __restrict__ PCrows, const int* __restrict__ PCcolumns, const Float* __restrict__ PCvalues, 
    Float inneriteration_tolerance,
    int inneriteration_print_modulo
);
 
#endif
