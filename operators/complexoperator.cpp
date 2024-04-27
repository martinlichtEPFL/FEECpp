
#include "complexoperator.hpp"

#include <cmath>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"
#include "simpleoperators.hpp"


FloatVector RealPart( const FloatVector& vec )
{
    assert( vec.getdimension() % 2 == 0 );
    return vec.getslice( 0, vec.getdimension() / 2 );
}

FloatVector ImaginaryPart( const FloatVector& vec )
{
    assert( vec.getdimension() % 2 == 0 );
    return vec.getslice( vec.getdimension() / 2, vec.getdimension() / 2 );
}

FloatVector ComplexFloatVector( const FloatVector& real, const FloatVector& imag )
{
    assert( real.getdimension() == imag.getdimension() );
    FloatVector ret( 2 * real.getdimension() );
    ret.setslice(                   0, real );
    ret.setslice( real.getdimension(), imag );
    return ret;
}





ComplexOperator::~ComplexOperator()
{
    if( part_real != nullptr && managing_real ) delete part_real;
    if( part_imag != nullptr && managing_imag ) delete part_imag;
}

void ComplexOperator::check() const  
{
    LinearOperator::check();    
    assert( part_real->getdimin()  == part_imag->getdimin()  );
    assert( part_real->getdimout() == part_imag->getdimout() );
}

std::string ComplexOperator::text() const 
{
    return text( false ); 
}

std::string ComplexOperator::text( const bool embellish ) const 
{
    return "Complex operator";
}

void ComplexOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    const auto src_real = RealPart( src );
    const auto src_imag = ImaginaryPart( src );

    const auto dest_real = scaling * ( *part_real * src_real - *part_imag * src_imag );
    const auto dest_imag = scaling * ( *part_real * src_imag + *part_imag * src_real );

    dest = ComplexFloatVector( dest_real, dest_imag );
}

