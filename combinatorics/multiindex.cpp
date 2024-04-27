
#include <sstream>
#include <string>
#include <vector>

#include "../basic.hpp"

#include "indexrange.hpp"

#include "multiindex.hpp"



MultiIndex::MultiIndex( const IndexRange& ir )
: IndexMap( ir, NonNegativeIntegers, std::vector<int>( ir.cardinality(), 0 ) )
{
    MultiIndex::check();
}

MultiIndex::MultiIndex( const IndexRange& ir, const std::vector<int>& vals )
: IndexMap( ir, NonNegativeIntegers, vals )
{
    assert( ir.cardinality() == vals.size() );
    MultiIndex::check();
}


MultiIndex::MultiIndex( const IndexRange& ir, const std::function<int(int)>& generator )
: IndexMap( ir, NonNegativeIntegers, generator )
{
    MultiIndex::check();
}


MultiIndex::MultiIndex( const IndexRange& ir, const std::initializer_list<int>& il )
: IndexMap( ir, NonNegativeIntegers, il )
{
    MultiIndex::check();
}
        





void MultiIndex::check() const
{
    #ifdef NDEBUG
    return;
    #endif
    
    IndexMap::check();
}

std::string MultiIndex::text( bool embellish ) const
{
    std::ostringstream ss;
    
    check();
    for( int p : getIndexRange() )
        ss << at( p ) << "\t";
    if( embellish ) 
        ss << nl;
    
    return ss.str();
}

// void MultiIndex::print( std::ostream& os, bool embellish ) const
// {
//     check();
//     for( int p : getIndexRange() )
//         os << at( p ) << "\t";
//     if( embellish ) 
//         os << nl;
// }

IndexRange MultiIndex::getIndexRange() const
{
    check();
    return IndexMap::getSourceRange();
}

const std::vector<int>& MultiIndex::getvalues() const
{
    check();
    return IndexMap::getvalues();
}
        





// const int& MultiIndex::at( int p ) const 
// {
//     check();
//     assert( getSourceRange().contains(p) );
//     return getvalues().at( getSourceRange().element2position(p) );
// }
// 
// int& MultiIndex::at( int p )
// {
//     check();
//     assert( getSourceRange().contains(p) );
//     return getvalues().at( getSourceRange().element2position(p) );
// }
// 
// const int& MultiIndex::operator[]( int p ) const 
// {
//     check();
//     assert( getSourceRange().contains(p) );
//     return getvalues().at( getSourceRange().element2position(p) );
// }
// 
// int& MultiIndex::operator[]( int p )
// {
//     check();
//     assert( getSourceRange().contains(p) );
//     return getvalues().at( getSourceRange().element2position(p) );
// }






int MultiIndex::absolute() const
{
    check();
    int ret = 0;
    for( int p : getIndexRange()  )
        ret += ::absolute<int>( at( p ) );
    return ret;
}

int MultiIndex::factorial() const
{
    check();
    
    assert( absolute() <= largest_factorial_base<int>() );
    assert( absolute() <= 20 );

    int ret = 1;
    for( int p : getIndexRange()  )
        ret *= factorial_integer( at( p ) );
    return ret;
}

Float MultiIndex::factorial_numerical() const
{
    check();
    
    Float ret = 1;
    for( int p : getIndexRange()  )
        ret *= ::factorial_numerical( at( p ) );
    return ret;
}

int MultiIndex::min() const
{
    check();
    
    int ret = getIndexRange().max();
    
    for( int p : getIndexRange()  )
        if( at(p) > 0 && p <= ret) 
            ret = p;
    
    assert( getIndexRange().contains(ret) );
    assert( at(ret) > 0 );
    for( int p : getIndexRange()  )
        assert( at(p) == 0 || p >= ret );
    
    return ret;
}

int MultiIndex::max() const
{
    check();
    
    int ret = getIndexRange().min();
    
    for( int p : getIndexRange()  )
        if( at(p) > 0 && p >= ret) 
            ret = p;
    
    assert( getIndexRange().contains(ret) );
    assert( at(ret) > 0 );
    for( int p : getIndexRange()  )
        assert( at(p) == 0 || p <= ret );
    
    return ret;
}






void MultiIndex::add( int p )
{
    check();
    assert( getIndexRange().contains(p) );
    at( p )++;
}

void MultiIndex::sub( int p )
{
    check();
    assert( getIndexRange().contains(p) );
    assert( at( p ) > 0 );
    at( p )--;
}



void MultiIndex::add( int p, int n )
{
    check();
    assert( getIndexRange().contains(p) );
    at( p ) += n;
}

void MultiIndex::sub( int p, int n )
{
    check();
    assert( getIndexRange().contains(p) );
    assert( at( p ) >= n );
    at( p ) -= n;
}




void MultiIndex::add( const MultiIndex& mi )
{
    check();
    mi.check();
    assert( getIndexRange() == mi.getIndexRange() );
    assert( is_comparable_with( mi ) );
    for( int p : getIndexRange() )
        add( p, mi.at( p ) );
    check();
}

void MultiIndex::sub( const MultiIndex& mi )
{
    check();
    mi.check();
    assert( getIndexRange() == mi.getIndexRange() );
    assert( is_comparable_with( mi ) );
    for( int p : getIndexRange() )
        sub( p , mi.at( p ) );
    check();
}






bool MultiIndex::is_comparable_with( const MultiIndex& mi ) const 
{
    check();
    mi.check();
    return getIndexRange() == mi.getIndexRange();
}


bool MultiIndex::is_less_than( const MultiIndex& mi ) const 
{
    check();
    mi.check();
    assert( is_comparable_with( mi ) );
    for( int p : getIndexRange() )
        if( at( p ) >= mi.at( p ) )
            return false;
    return true;
}

bool MultiIndex::is_equal_to( const MultiIndex& mi ) const
{
    check();
    mi.check();
    assert( is_comparable_with( mi ) );
    for( int p : getIndexRange() )
        if( at( p ) != mi.at( p ) )
            return false;
    return true;
}

