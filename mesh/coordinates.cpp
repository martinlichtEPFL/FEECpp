

// #include <cassert>
#include <algorithm>
#include <istream>
#include <ostream>
#include <sstream>
#include <iterator>
#include <vector>

#include "coordinates.hpp"


Coordinates::Coordinates( int dimension, int number )
: 
    dimension( dimension ), 
    number( number ), 
    data( dimension * number )
{
    Coordinates::check();
}

Coordinates::Coordinates( int dimension, int number, const std::vector<Float>& data )
: 
    dimension( dimension ), 
    number( number ), 
    data( data )
{
    Coordinates::check();
}

Coordinates::~Coordinates()
{
    // Coordinates::check();
}

void Coordinates::check() const
{
    assert( dimension >= 0 && number >= 0 );
    if( data.size() != 0 )
        assert( data.size() == dimension * number );
}



// void Coordinates::print( std::ostream& os ) const 
// {
//     os << text();
// }

std::string Coordinates::text() const
{
    std::stringstream os;
    
    os << "dimension: " << dimension << " - #vertices: " << number << nl;
    for( int n = 0; n < number; n++ ) {
        for( int d = 0; d < dimension; d++ )
            os << getdata( n, d ) << " ";
        os << nl;
    }

    return os.str();
}


// void Coordinates::read( std::istream& is ) 
// {
//     for( int n = 0; n < number; n++ ) {
//         for( int d = 0; d < dimension; d++ ) {
//             Float temp;
//             is >> temp;
//             setdata( n, d, temp );
//             assert( ! is.fail() );
//         }
//     }
// }





int Coordinates::getdimension() const
{
    return dimension;
}

int Coordinates::getnumber() const
{
    return number;
}

IndexRange Coordinates::getIndexRange() const
{
    return IndexRange( 0, getnumber() - 1 );
}


Float Coordinates::getdata( int n, int d ) const
{
    assert( 0 <= n && n < number && 0 <= d && d < dimension );
    return data.at( n * dimension + d );
}

void Coordinates::setdata( int n, int d, Float v )
{
    assert( 0 <= n && n < number && 0 <= d && d < dimension );
    data.at( n * dimension + d ) = v;
}


FloatVector Coordinates::getvectorclone( int n ) const
{
    assert( 0 <= n && n < number );
    return getvectorclone( n, 1. );
}

/* get range of coordinates */
        
Float Coordinates::getmin( int d ) const {
    assert( 0 <= d && d < dimension );
    Float ret = getdata(0,d);
    for( int n = 1; n < number; n++ )
        ret = minimum(ret,getdata(n,d));
    return ret;
}
Float Coordinates::getmax( int d ) const {
    assert( 0 <= d && d < dimension );
    Float ret = getdata(0,d);
    for( int n = 1; n < number; n++ )
        ret = maximum(ret,getdata(n,d));
    return ret;
}
        
        
        
FloatVector Coordinates::getvectorclone( int n, Float scale ) const
{
    assert( 0 <= n && n < number );
    FloatVector ret( dimension );
    for( int d = 0; d < dimension; d++ )
        ret[ d ] = scale * data.at( n * dimension + d );
    return ret;
}

void Coordinates::loadvector( int n, const FloatVector& input ) 
{
    loadvector( n, input, 1. );
}

void Coordinates::loadvector( int n, const FloatVector& input, Float scale ) 
{
    assert( 0 <= n && n < number );
    assert( input.getdimension() == dimension );
    for( int d = 0; d < dimension; d++ )
        data.at( n * dimension + d ) = scale * input[ d ];
}



/* get/set coordinates as vectors  */
        
FloatVector Coordinates::getdimensionclone( int d, Float s ) const 
{
    assert( 0 <= d && d < dimension );
    FloatVector ret( number );
    for( int n = 0; n < number; n++ )
        ret[ n ] = s * data.at( n * dimension + d );
    return ret;    
}

void Coordinates::loaddimension( int d, const FloatVector& value, Float s )
{
    assert( 0 <= d && d < dimension );
    assert( value.getdimension() == number );
    for( int n = 0; n < number; n++ )
        data.at( n * dimension + d ) = s * value[ n ];    
}





void Coordinates::scale( Float alpha )
{
    for( int n = 0; n < number; n++ )
        for( int d = 0; d < dimension; d++ )
            data.at( n * dimension + d ) *= alpha;
}
                                
void Coordinates::scale( FloatVector alphas )
{
    for( int n = 0; n < number; n++ )
        for( int d = 0; d < dimension; d++ )
            data.at( n * dimension + d ) *= alphas[d];
}
                                
void Coordinates::shift( const FloatVector& add )
{
    assert( add.getdimension() == dimension );
    for( int n = 0; n < number; n++ ) {
        FloatVector temp = getvectorclone( n );
        temp += add;
        loadvector( n, temp );
    }
}

void Coordinates::lineartransform( const LinearOperator& op )
{
    assert( op.getdimin() == dimension );
    assert( op.getdimout() == dimension );
    for( int n = 0; n < number; n++ ) {
        FloatVector temp = getvectorclone( n );
        temp = op * temp;
        loadvector( n, temp );
    }
}



void Coordinates::append( const Coordinates& co )
{
    assert( dimension == co.getdimension() );
    data.insert( this->data.end(), co.data.begin(), co.data.end() );
    number += co.number;
}
                
void Coordinates::append( const FloatVector& v )
{
    assert( dimension == v.getdimension() );
    const std::vector<Float>& t = v.getdata();
    data.insert( this->data.end(), t.begin(), t.end() );
    number++;
}


void Coordinates::addcapacity( int additional_capacity )
{
    assert( additional_capacity >= 0 );
    data.reserve( data.size() + dimension * additional_capacity );
}

void Coordinates::addcoordinates( int add_number )
{
    assert( add_number >= 0 );
    data.resize( data.size() + dimension * add_number );
    number += add_number;
}





DenseMatrix Coordinates::getLinearPart( const IndexMap& im ) const
{
    assert( im.getTargetRange() == getIndexRange() );
    
    IndexRange imsrc = im.getSourceRange();
    assert( imsrc.min() == 0 && imsrc.max() <= getdimension() );
    
    DenseMatrix ret( getdimension(), maximum(0,imsrc.max()-1) );
    assert( ret.getdimout() == getdimension() );
    assert( ret.getdimin() == maximum(0,imsrc.max()-1) );
    
    for( int p = 1; p <= imsrc.max(); p++ )
        ret.setcolumn( p-1, 
            getvectorclone( im[p] ) - getvectorclone( im[0] ) 
            );
    
    return ret;
}

FloatVector Coordinates::getShiftPart( const IndexMap& im ) const
{
    assert( im.getTargetRange() == getIndexRange() );
    IndexRange imsrc = im.getSourceRange();
    assert( !(im.isempty()) && imsrc.min() == 0 && imsrc.max() <= getdimension() );
    int index = im[0];
    return getvectorclone( index );
}




FloatVector Coordinates::getCenter() const
{
    FloatVector center( dimension, 0. );
    for( int i = 0; i < number * dimension; i++ )
        center[ i % dimension ] += data[ ( i / dimension ) * dimension + ( i % dimension ) ];
    for( int d = 0; d < dimension; d++ )
        center[ d ] /= number;
    return center;
}




std::vector<Float>& Coordinates::raw()
{
    return data;
}

const std::vector<Float>& Coordinates::raw() const
{
    return data;
}




        


std::size_t Coordinates::memorysize() const
{
    return sizeof(*this) + data.size() * sizeof(Float);
}




bool Coordinates::is_equal_to( const Coordinates& coords_left, const Coordinates& coords_right )
{
    if( coords_left.getnumber() != coords_right.getnumber() )
      return false;
    
    if( coords_left.getdimension() != coords_right.getdimension() )
      return false;
    
    for( int n = 0; n < coords_left.getnumber(); n++ )
    for( int d = 0; d < coords_left.getdimension(); d++ )
      if( not is_numerically_close( coords_left.getdata( n, d ), coords_right.getdata( n, d ), 0.1 ) )
        return false;
    
    return true;
}













