
#include "simpleoperators.hpp"

#include <cmath>

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"


IdentityOperator::IdentityOperator( int dimension )
: LinearOperator( dimension )
{
    IdentityOperator::check();
}

IdentityOperator::~IdentityOperator()
{
    IdentityOperator::check();
}

void IdentityOperator::check() const  
{
    LinearOperator::check();
}

std::string IdentityOperator::text() const 
{
    return "Identity Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin());
}

void IdentityOperator::apply( FloatVector& dest, const FloatVector& src, Float s ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    for( int p = 0; p < getdimin(); p++ )
        dest.setentry( p, s * src.getentry( p ) );
        
}














ZeroOperator::ZeroOperator( int dimension )
: LinearOperator( dimension, dimension )
{
    ZeroOperator::check();
}

ZeroOperator::ZeroOperator( int dimout, int dimin )
: LinearOperator( dimout, dimin )
{
    ZeroOperator::check();
}

ZeroOperator::~ZeroOperator()
{
    ZeroOperator::check();
}

void ZeroOperator::check() const  
{
    LinearOperator::check();
}

std::string ZeroOperator::text() const 
{
    return "Zero Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin());
}

void ZeroOperator::apply( FloatVector& dest, const FloatVector& src, Float s ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    dest.clear();
}














ScalingOperator::ScalingOperator( int dimension, Float s )
: LinearOperator( dimension ), 
  scaling(s)
{
    ScalingOperator::check();
}

ScalingOperator::~ScalingOperator()
{
    ScalingOperator::check();
}

void ScalingOperator::check() const  
{
    LinearOperator::check();
}

std::string ScalingOperator::text() const 
{
    return "Scaling Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + ": s = " + std::to_string( scaling );
}


void ScalingOperator::setscaling( Float s )
{
    scaling = s;
}

Float ScalingOperator::getscaling() const
{
    return scaling;
}

void ScalingOperator::apply( FloatVector& dest, const FloatVector& src, Float s ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    for( int p = 0; p < getdimin(); p++ )
        dest.setentry( p, s * scaling * src.getentry( p ) );
        
}














DiagonalOperator::DiagonalOperator( int dimension, Float scale )
: LinearOperator( dimension ), 
  diagonal( FloatVector( dimension, scale ) )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( const FloatVector& dia )
: LinearOperator( dia.getdimension() ), 
  diagonal( dia )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( FloatVector&& dia )
: LinearOperator( dia.getdimension() ), 
  diagonal( std::move(dia) )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( int dimension, const ScalingOperator& scaling )
: LinearOperator( dimension ), 
  diagonal( FloatVector( dimension, scaling.getscaling() ) )
{
    DiagonalOperator::check();
}

DiagonalOperator::DiagonalOperator( int dimension, const std::function<Float(int)>& generator )
: LinearOperator( dimension ), 
  diagonal( FloatVector( dimension, generator ) )
{
    DiagonalOperator::check();
}
        

DiagonalOperator::~DiagonalOperator()
{
    /* Nothing */ 
}

void DiagonalOperator::check() const  
{
    LinearOperator::check();    
    diagonal.check();
    assert( getdimin() == getdimout() );
    assert( getdimin() == diagonal.getdimension() );
}

std::string DiagonalOperator::text() const 
{
    return text( false ); 
}

std::string DiagonalOperator::text( const bool embellish ) const 
{
    return "Diagonal operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) 
            + ( embellish ? "\n" + tab_each_line( diagonal.text() ) : "" );
}


FloatVector& DiagonalOperator::getdiagonal()
{
    return diagonal;
}

const FloatVector& DiagonalOperator::getdiagonal() const
{
    return diagonal;
}

void DiagonalOperator::apply( FloatVector& dest, const FloatVector& src, Float scaling ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );
    
    for( int p = 0; p < getdimout(); p++ ) {
        dest.setentry( p, scaling * diagonal.at(p) * src.getentry( p ) );
    }
    
}

const DiagonalOperator DiagonalOperator::sqrt() const
{
    return DiagonalOperator( getdimin(), [this](int i) -> Float {
        assert( this->diagonal[i] >= 0 ); 
        return std::sqrt( diagonal[i] );
    } );
}















LambdaOperator::LambdaOperator( int n, std::function<FloatVector(const FloatVector&)> func )
: LambdaOperator( n, n, func )
{
    LambdaOperator::check();
}

LambdaOperator::LambdaOperator( int dimout, int dimin, std::function<FloatVector(const FloatVector&)> func )
: LinearOperator( dimout, dimin ), func(func)
{
    LambdaOperator::check();
}

LambdaOperator::~LambdaOperator()
{
    LambdaOperator::check();
}

void LambdaOperator::check() const  
{
    LinearOperator::check();
}

std::string LambdaOperator::text() const 
{
    return "Lambda Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin());
}

void LambdaOperator::apply( FloatVector& dest, const FloatVector& src, Float s ) const 
{
    check();
    src.check();
    dest.check();
    
    assert( getdimin() == getdimout() );
    assert( getdimin() == src.getdimension() );
    assert( getdimout() == dest.getdimension() );

    dest = s * func( src );
        
}
