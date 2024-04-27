#ifndef INCLUDEGUARD_SPARSE_CSRRAINBOW_HPP
#define INCLUDEGUARD_SPARSE_CSRRAINBOW_HPP

#include <string>
#include <vector>
#include <utility>

#include "../basic.hpp"
#include "matcsr.hpp"



struct Rainbow
{
    explicit Rainbow( const MatrixCSR&, bool do_shuffle = false );

    void check() const;

    std::string text() const;

    // number of colors used 
    int num_rows;
    int num_colors;

    // for each row index, give the color index 
    std::vector<int> F; 

    // CSR-like construction 
    std::vector<int> B; 
    std::vector<int> R;
};


#endif
