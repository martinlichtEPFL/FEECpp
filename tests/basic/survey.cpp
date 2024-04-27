
#include <cstdio>
#include <cfloat>
#include "../../basic.hpp"

using namespace std;

int main( int argc, char *argv[] )
{
    printf("Survey of machine data (integral and floating-point)");

    /* output integer parameters */

    printf("    size of char      : %ju\n", (uintmax_t)sizeof(       char) );
    printf("    size of short int : %ju\n", (uintmax_t)sizeof(  short int) );
    printf("    size of int       : %ju\n", (uintmax_t)sizeof(        int) );
    printf("    size of long int  : %ju\n", (uintmax_t)sizeof(   long int) );
    printf("    size of long long : %ju\n", (uintmax_t)sizeof(  long long) );
    
    printf("    size of size_t    : %ju\n", (uintmax_t)sizeof(     size_t) );
    printf("    size of size_t    : %ju\n", (uintmax_t)sizeof(std::size_t) );
    printf("    SIZE_MAX          : %ju\n", (uintmax_t) SIZE_MAX           );

    /* output floating-point parameters */

    printf("    \nsizes of floating-point types\n");
    printf("    size of float       : %ju\n", (uintmax_t)sizeof(      float) );
    printf("    size of double      : %ju\n", (uintmax_t)sizeof(     double) );
    printf("    size of long double : %ju\n", (uintmax_t)sizeof(long double) );
    
    printf("    \nFloating-point mode flags\n");
    printf("    FLT_ROUNDS      : %ld\n", (long) FLT_ROUNDS      );
    printf("    FLT_EVAL_METHOD : %ld\n", (long) FLT_EVAL_METHOD );
    printf("    FLT_RADIX       : %ld\n", (long) FLT_RADIX       );
    printf("    DECIMAL_DIG     : %ld\n", (long) DECIMAL_DIG     );

    printf("    \nFloating-point properties\n");
    #if __cplusplus >= 201703L
    printf("    FLT/DBL/LDBL_DECIMAL_DIG : %ld\t%ld\t%ld\n", (long) FLT_DECIMAL_DIG, (long) DBL_DECIMAL_DIG, (long) LDBL_DECIMAL_DIG );
    #endif //__cplusplus >= 201703L 
    printf("    FLT/DBL/LDBL_DIG         : %ld\t%ld\t%ld\n", (long) FLT_DIG,         (long) DBL_DIG,         (long) LDBL_DIG         );
    printf("    FLT/DBL/LDBL_MANT_DIG    : %ld\t%ld\t%ld\n", (long) FLT_MANT_DIG,    (long) DBL_MANT_DIG,    (long) LDBL_MANT_DIG    );
    printf("    FLT/DBL/LDBL_MIN_EXP     : %ld\t%ld\t%ld\n", (long) FLT_MIN_EXP,     (long) DBL_MIN_EXP,     (long) LDBL_MIN_EXP     );
    printf("    FLT/DBL/LDBL_MAX_EXP     : %ld\t%ld\t%ld\n", (long) FLT_MAX_EXP,     (long) DBL_MAX_EXP,     (long) LDBL_MAX_EXP     );
    printf("    FLT/DBL/LDBL_MIN_10_EXP  : %ld\t%ld\t%ld\n", (long) FLT_MIN_10_EXP,  (long) DBL_MIN_10_EXP,  (long) LDBL_MIN_10_EXP  );
    printf("    FLT/DBL/LDBL_MAX_10_EXP  : %ld\t%ld\t%ld\n", (long) FLT_MAX_10_EXP,  (long) DBL_MAX_10_EXP,  (long) LDBL_MAX_10_EXP  );
    #if __cplusplus >= 201703L 
    printf("    FLT/DBL/LDBL_HAS_SUBNORM : %ld\t%ld\t%ld\n", (long) FLT_HAS_SUBNORM, (long) DBL_HAS_SUBNORM, (long) LDBL_HAS_SUBNORM );
    #endif //__cplusplus >= 201703L 

    printf("    \nFloating-point Minima and Maxima\n");
    printf("    FLT/DBL/LDBL_MIN      : %Le\t%Le\t%Le\n", (long double) FLT_MIN,      (long double) DBL_MIN,      (long double) LDBL_MIN      );
    printf("    FLT/DBL/LDBL_MAX      : %Le\t%Le\t%Le\n", (long double) FLT_MAX,      (long double) DBL_MAX,      (long double) LDBL_MAX      );
    printf("    FLT/DBL/LDBL_EPSILON  : %Le\t%Le\t%Le\n", (long double) FLT_EPSILON,  (long double) DBL_EPSILON,  (long double) LDBL_EPSILON  );
    #if __cplusplus >= 201703L 
    printf("    FLT/DBL/LDBL_TRUE_MIN : %Le\t%Le\t%Le\n", (long double) FLT_TRUE_MIN, (long double) DBL_TRUE_MIN, (long double) LDBL_TRUE_MIN );
    #endif //__cplusplus >= 201703L 


    printf("\nProject-defined floating-point type\n");
    printf("Float size:             %ju\n", (uintmax_t)sizeof(Float) );
    
    printf("Float decimal digits:   %ld\n",  (long int) std::numeric_limits<Float>::digits10         );
    printf("Float digits:           %ld\n",  (long int) std::numeric_limits<Float>::digits           );
    printf("Float min exponent:     %ld\n",  (long int) std::numeric_limits<Float>::min_exponent     );
    printf("Float max exponent:     %ld\n",  (long int) std::numeric_limits<Float>::max_exponent     );
    printf("Float min exponent 10:  %ld\n",  (long int) std::numeric_limits<Float>::min_exponent10   );
    printf("Float max exponent 10:  %ld\n",  (long int) std::numeric_limits<Float>::max_exponent10   );
    printf("Float has denorm:       %ld\n",  (long int) std::numeric_limits<Float>::has_denorm       );
    printf("Float has denorm loss:  %ld\n",  (long int) std::numeric_limits<Float>::has_denorm_loss  );

    printf("Float denormalized min: %Le\n",  (long double) std::numeric_limits<Float>::denorm_min()  );
    printf("Float minimum:          %Le\n",  (long double) std::numeric_limits<Float>::min()         );
    printf("Float maximum:          %Le\n",  (long double) std::numeric_limits<Float>::max()         );
    printf("Float machine epsilon:  %Le\n",  (long double) std::numeric_limits<Float>::epsilon()     );

    printf("Float rounding error:    %Lf\n", (long double) std::numeric_limits<Float>::round_error() );
    printf("Float has quiet NaN:     %d\n",  (int) std::numeric_limits<Float>::has_quiet_NaN         );
    printf("Float has signaling NaN: %d\n",  (int) std::numeric_limits<Float>::has_signaling_NaN     );
    printf("Float is IEC-559:        %d\n",  (int) std::numeric_limits<Float>::is_iec559             );
    printf("Float detectes tinyness: %d\n",  (int) std::numeric_limits<Float>::tinyness_before       );
    
    printf("Machine epsilon variable:         %.25Le\n", (long double)machine_epsilon );
    printf("Machine epsilon variable (sqrt ): %.25Le\n", (long double)std::sqrt(machine_epsilon) );
    printf("Machine epsilon variable (Sqrt ): %.25Le\n", (long double)Sqrt(machine_epsilon) );
    // printf("Machine epsilon variable (Sqrt_): %.25Le\n", (long double)Sqrt_(machine_epsilon) );
    printf("Desired precision variable:      %Le\n",     (long double)desired_precision );
    printf("Desired closeness variable:      %Le\n",     (long double)desired_closeness );

    // printf("test: %.25Le\n", (long double)std::sqrt(4 * machine_epsilon * machine_epsilon) / 2 );
    // printf("test: %.25Le\n", (long double)std::sqrt(machine_epsilon * machine_epsilon)/machine_epsilon );
    


    // {
    //     unsigned int i = 1000;
    //     long double a = std::numeric_limits<long double>::epsilon();
    //     long double x = a;
    //     while ( i --> 0 ) { 
    //         x = ( x + a / x ) / 2.;
    //         printf("i=%u x=%Le\n", i, (long double)x );
    //     }
    // }
    
    assert( machine_epsilon < 1e-10 );

    

    return 0;
}
