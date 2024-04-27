#ifndef INCLUDEGUARD_UTILITY_RANDOM_HPP
#define INCLUDEGUARD_UTILITY_RANDOM_HPP


#include "../basic.hpp"

void seed_random_integer();

unsigned int random_integer();

unsigned int flip_coin( Float prob_zero = 0.5 );

Float random_uniform();

Float gaussrand();



#endif
