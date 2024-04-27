

#include "linearoperator.hpp"

#include "floatvector.hpp"

LinearOperator::LinearOperator( int out, int in )
: dimout( out ), dimin( in )
{
    LinearOperator::check();
}

LinearOperator::LinearOperator( int dim )
: dimout( dim ), dimin( dim )
{
    LinearOperator::check();
}

LinearOperator::~LinearOperator()
{
    LinearOperator::check();
}

int LinearOperator::getdimout() const
{
    /* No check here */
    return dimout;
}

int LinearOperator::getdimin() const
{
    /* No check here */
    return dimin;
}

void LinearOperator::check() const
{
    assert( dimout >= 0 && dimin >= 0 );
}

// void LinearOperator::print( std::ostream& os ) const
// {
//     os << text();
// }

bool LinearOperator::issquare() const
{
    check();
    return getdimin() == getdimout();
}

/* x := A y */
FloatVector LinearOperator::apply( const FloatVector& src, Float scaling ) const
{
    check();
    FloatVector ret = createoutputvector();
    apply( ret, src, scaling );
    return ret;
}

FloatVector LinearOperator::createinputvector() const
{
    return FloatVector( getdimin() );
}

FloatVector LinearOperator::createoutputvector() const
{
    return FloatVector( getdimout() );
}





