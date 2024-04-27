
#include <algorithm>
#include <fstream>
#include <istream>
#include <map>
#include <ostream>
#include <string>
#include <vector>
#include <utility>



#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "io.coordinates.hpp"




void writeCoordinates( const char* filename, const Coordinates& coords, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeCoordinates( myfile, coords, sugar );
    myfile.close();
}

Coordinates readCoordinates( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    Coordinates coords = readCoordinates( myfile );
    myfile.close();
    return coords;
}


void writeCoordinates( std::ostream& out, const Coordinates& coords, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Writing coordinates" << nl;
    if( sugar ) out << "Dimension: "         << nl;
    out << coords.getdimension() << nl;
    if( sugar ) out << "Number of points: "  << nl;
    out << coords.getnumber() << nl;
    
    /* data */
    for( int p = 0; p < coords.getnumber(); p++ ) {
        if( sugar ) out << p << ": ";
        for( int d = 0; d < coords.getdimension(); d++ )
            out << coords.getdata(p,d) << " ";
        out << nl;
    }
    
    assert( out.good() );
    
}

Coordinates readCoordinates( std::istream& in )
{
    int dimension;
    int numpoints;
    in >> dimension >> numpoints;
    Coordinates coords( dimension, numpoints );
    for( int p = 0; p < numpoints; p++ ) {
        for( int d = 0; d < dimension; d++ ) {
            Float tempvalue;
            in >> tempvalue;
            coords.setdata( p, d, tempvalue );
        }
    }
    
    assert( in.good() );
    
    return coords;
}



        
        
        
        
        
