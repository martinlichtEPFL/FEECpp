// g++ sorthacktest.cpp -o sorthacktest

#include <cassert>
#include <cstdlib>

#include <algorithm>
#include <vector>

#include "../../basic.hpp"



int main( int argc, char *argv[] )
{
    
    LOG << "This is a line" << nl;
    NOTE "This is a note";
    PING;
    PING;
    PING;
    LOG << "This is a test! " << 5 << nl;
    LOGPRINTF( "%i%c%d\n", 1, '-', 3 );
    LOGPRINTF( "0.123456789e-7: %e\n",   (double)     0.1234567890123456789e-7 );
    LOGPRINTF( "0.123456789e-7: % Le\n", (long double)0.1234567890123456789e-7 );
    PING;
    // openContext();
    NOTE "";
    NOTE "Notice";
    WARNING "Warning";
    PING;
    // closeContext();
    return 0;
}
