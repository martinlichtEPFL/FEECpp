#include <cstdio>
#include <cmath>
#include <cfloat>

#include <limits>

#include "../../basic.hpp"

typedef struct { float upper; float lower; } worst_and_worst;

static worst_and_worst calculate_sqrt_distance( float x ) {
    // Example function that just prints the float
    float s1 = std::sqrt(x);
    float s2 = Sqrt(x,200);
    return { s2/s1, s1/s2 };
    // std::printf( "%c %.30e %.10e %.10e %.20e\n", std::isnormal(x) ? 'n' : 's', x,   s2,   s1, s2/s1 );
}

int main() 
{

    LOGPRINTF( "%e %e %e \n", FLT_MIN, FLT_TRUE_MIN, FLT_MIN / 2. );



    static_assert( sizeof(float) == sizeof(int) );

    union {
        float f;
        unsigned int i;
    } u;

    unsigned int maxUint = 10000; // std::numeric_limits<unsigned int>::max();

    u.f = std::numeric_limits<float>::epsilon();
    calculate_sqrt_distance( u.f );
    
    // return 0;

    float worst_upper = 0.;
    float worst_lower = 0.;
    
    #if defined(_OPENMP)
    #pragma omp parallel for reduction(max:worst_upper,worst_lower)
    #endif
    for( int i = 0; i < maxUint; ++i )
    {
        u.i = i;

        // Optionally, break or continue on specific conditions to avoid NaNs, infinities, etc.
        if (std::isnan(u.f) || std::isinf(u.f)) {
            continue; // Skip this iteration
        }

        if( u.f < 0 ) 
            continue; // skip from here on 

        auto worst = calculate_sqrt_distance(u.f);

        worst_upper = std::max( worst_upper, worst.upper );
        worst_lower = std::max( worst_lower, worst.lower );

        // LOGPRINTF( "%.10e %.10e %.10e %.10e \n", worst_upper, worst.upper, worst_lower, worst.lower );
    
    }

    // std::printf( "%c %.30e %.10e %.10e %.20e\n", std::isnormal(x) ? 'n' : 's', x,   s2,   s1, s2/s1 );
    
    // Make sure to also call foo() on the max float value
    u.i = maxUint;
    calculate_sqrt_distance(u.f);

    LOGPRINTF( "%e %e \n", worst_upper, worst_lower );

    return 0;
}
