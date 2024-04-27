
#include <algorithm>
#include <sstream>

#include "indexrange.hpp"
#include "indexmap.hpp"





// IndexMap::IndexMap( const IndexRange& range )
// : IndexMap( range, range )
// {
//     check();
// }
// 
// IndexMap::IndexMap( const IndexRange& from, const IndexRange& to )
// : src(from), dest(to), values( maximum( src.max() - src.min() + 1, 0 ), to.min() )
// {
//     if( src.max() >= src.min() )
//       LOG << "Index Map initialized without actual values" << nl;
//     src.check();
//     to.check();
//     check();
// }

IndexMap::IndexMap( const IndexRange& range, const std::vector<int>& values )
: IndexMap( range, range, values )
{}

IndexMap::IndexMap( const IndexRange& from, const IndexRange& to, const std::vector<int>& values )
: src(from), dest(to), values(values)
{}

IndexMap::IndexMap( const IndexRange& range, const std::function<int(int)>& generator )
: IndexMap( range, range, generator )
{}

IndexMap::IndexMap( const IndexRange& from, const IndexRange& to, const std::function<int(int)>& generator )
: src(from), dest(to), values()
{
    if( from.isempty() ) {
    
        values.resize(0);
    
    } else {
    
        values.reserve( maximum( src.max() - src.min() + 1, 0 ) );
        for( int e = src.min(); e <= src.max(); e++ )
            // values.at( src.element2position(e) ) = generator(e);
            values.emplace_back( generator(e) );
        
    }
    IndexMap::check();
}

IndexMap::IndexMap( const IndexRange& range, const std::initializer_list<int>& values )
: IndexMap( range, std::vector<int>(values) )
{}

IndexMap::IndexMap( const IndexRange& from, const IndexRange& to, const std::initializer_list<int>& values )
: IndexMap( from, to, std::vector<int>(values) )
{}

IndexMap::IndexMap( const IndexRange& from, const IndexRange& to, const int& value )
: IndexMap( from, to, [=](int) -> int { return value; } )
{}


        
        



void IndexMap::check() const 
{ 
    #ifdef NDEBUG
    return;
    #endif
    
    src.check();
    dest.check();
    
    assert( getSourceRange().cardinality() == values.size() );
    
    if( values.size() > 0 ) {
        
        assert( ! getSourceRange().isempty() );
        
        assert( maximum( src.max() - src.min() + 1, 0 ) == values.size() );
    
        assert( ! getTargetRange().isempty() );
        
        for( int a = src.min(); a <= src.max(); a++ )
            assert( dest.contains( values.at( a - src.min() ) ) );
        
    } else {
        
        assert( getSourceRange().isempty() );
        
    }
        
}

std::string IndexMap::text( bool embellish ) const 
{
    std::ostringstream ss;
    
    if( embellish ) {
        ss << ' ' << getSourceRange().text() << " => " << getTargetRange().text() << nl;
        for( int i : getSourceRange() )
            ss << '\t' << i << " -> " << at( i ) << nl;
    } else {
        ss << getSourceRange() << "\n" << getTargetRange() << nl;
        for( int i : getSourceRange() )
            ss << " " << at( i );
        ss << nl;
    }
    
    return ss.str();
}

// void IndexMap::print( std::ostream& os, bool embellish ) const 
// {
//     os << text( embellish );
// }

const IndexRange& IndexMap::getSourceRange() const 
{
    return src;
}

const IndexRange& IndexMap::getTargetRange() const 
{
    return dest;
}





int& IndexMap::at( int i ) &
{
    check();
    assert( !src.isempty() );
    assert( src.contains(i) );
    assert( 0 <= i - src.min() );
    assert( i - src.min() < values.size() );
    return values.at( i - src.min() );
}

const int& IndexMap::at( int i ) const&
{
    check();
    assert( !src.isempty() );
    assert( src.contains(i) );
    assert( 0 <= i - src.min() );
    assert( i - src.min() < values.size() );
    return values.at( i - src.min() );
}

int& IndexMap::operator[]( int i ) &
{
    check();
    assert( !src.isempty() );
    assert( src.contains(i) );
    assert( 0 <= i - src.min() );
    assert( i - src.min() < values.size() );
    return values[ i - src.min() ];
}

const int& IndexMap::operator[]( int i ) const&
{
    check();
    assert( !src.isempty() );
    assert( src.contains(i) );
    assert( 0 <= i - src.min() );
    assert( i - src.min() < values.size() );
    return values[ i - src.min() ];
}

const std::vector<int>& IndexMap::getvalues() const &
{
    return values;
}





bool IndexMap::isempty() const 
{
    check();
    return getSourceRange().isempty();
}


bool IndexMap::isinjective() const 
{
    check();
    
    if( getSourceRange().isempty() )
        return true;
    
    for( int a = src.min(); a <= src.max(); a++ )
        for( int b = src.min(); b <= src.max(); b++ )
            if( a != b && values.at( a - src.min() ) == values.at( b - src.min() ) )
                return false;
    
    return true;
}

bool IndexMap::issurjective() const 
{
    check();
    
    if( getSourceRange().isempty() )
        return true;
    
    for( int a = dest.min(); a <= dest.max(); a++ ) 
    {
        bool flag = false;
        for( int b = src.min(); b <= src.max(); b++ )
            if( values.at( b - src.min() ) == a )
                flag = true;
        if( !flag )
            return false;
    }
    
    return true;
}

bool IndexMap::isbijective() const 
{
    check();
    return isinjective() && issurjective();
}

bool IndexMap::isstrictlyascending() const
{
    check();
    
    if( getSourceRange().isempty() )
        return true;
    
    for( int a = src.min(); a < src.max(); a++ )
        if( values.at( a - src.min() ) >= values.at( a - src.min() + 1 ) )
            return false;
    
    return true;
}

bool IndexMap::has_value_in_range( int p ) const
{
    check();
    assert( getTargetRange().contains(p) );
    for( int i : src )
        if( at(i) == p )
            return true;
    return false;
} 

int IndexMap::preimageof( int p ) const
{
    check();
    assert( getTargetRange().contains(p) );
    for( int i : src )
        if( at(i) == p )
            return i;
    unreachable();
} 
        
bool IndexMap::is_comparable_with( const IndexMap& im ) const
{
    check();
    return getSourceRange() == im.getSourceRange() && getTargetRange() == im.getTargetRange();
} 

bool IndexMap::is_equal_to( const IndexMap& im ) const
{
    check();
    im.check();
    assert( is_comparable_with( im ) );
    return getvalues() == im.getvalues();
} 

bool IndexMap::is_less_than( const IndexMap& im ) const
{
    check();
    im.check();
    assert( is_comparable_with( im ) );
    return getvalues() < im.getvalues();
} 










IndexMap mergeSigmas( const IndexMap& left, const IndexMap& right, int& sign )
{
    sign = 1;
    if( left.getSourceRange().isempty()  ) return right;
    if( right.getSourceRange().isempty() ) return left;

    assert(  left.getSourceRange().min() == 1 );
    assert( right.getSourceRange().min() == 1 );
    assert(  left.getTargetRange().min() == 0 );
    assert( right.getTargetRange().min() == 0 );

    const int k =  left.getSourceRange().max();
    const int l = right.getSourceRange().max();

    const int n = maximum( left.getTargetRange().max(), right.getTargetRange().max() );

    std::vector<int> values(k+l);
    for( int i = 0; i < k; i++ ) values[  i] =  left.getvalues()[i];
    for( int i = 0; i < l; i++ ) values[k+i] = right.getvalues()[i];

    sign = 1;
    // lazy bubble sort, easy to count the swaps
    for( int i = 0; i < k+l; i++ )
    for( int j = 0; j < k+l; j++ )
    {
        if( i < j and values[i] > values[j] )
        {
            std::swap( values[i], values[j] );
            sign = -sign;
        }
    }

    for( int i = 1; i < k+l; i++ )
        if( values[i-1] == values[i] )
            sign = 0;

    return IndexMap( IndexRange(1,k+l), IndexRange(0,n), values );
}









IndexMap expand_zero( const IndexMap& im, int p )
{
    const auto& src_range = im.getSourceRange();
    const auto& dst_range = im.getTargetRange();
    
    assert( not dst_range.isempty() );
    assert( dst_range.min() == 0    );
    
    if( src_range.isempty() ) {
        
        return IndexMap( IndexRange(0,0), dst_range, {p} );
        
    } else {
        
        assert( src_range.min() == 0 );
        auto new_values = im.getvalues();
        new_values.push_back( p );
        std::sort( new_values.begin(), new_values.end() );
        const auto ret = IndexMap( IndexRange(0,src_range.max()+1), dst_range, new_values );
        assert( ret.isstrictlyascending() );
        return ret;
        
    }
}


IndexMap expand_one( const IndexMap& im, int p )
{
    const auto& src_range = im.getSourceRange();
    const auto& dst_range = im.getTargetRange();
    
    assert( not dst_range.isempty() );
    assert( dst_range.min() == 0    );
    
    if( src_range.isempty() ) {
        
        return IndexMap( IndexRange(1,1), dst_range, {p} );
        
    } else {
        
        assert( src_range.min() == 1 );
        auto new_values = im.getvalues();
        new_values.push_back( p );
        std::sort( new_values.begin(), new_values.end() );
        const auto ret = IndexMap( IndexRange(1,src_range.max()+1), dst_range, new_values );
        assert( ret.isstrictlyascending() );
        return ret;
        
    }
}



IndexMap complement_sigma( const IndexMap& sigma, const int n )
{
    const auto source_range = sigma.getSourceRange();
    const auto target_range = sigma.getTargetRange();

    assert( n >= 0 );
    assert( not target_range.isempty() );
    assert( target_range.min() == 0 );
    assert( target_range.max() == n );

    if( source_range.isempty() )
        return IndexMap( IndexRange(0,n), IndexRange(0,n), [](int i)->int{return i;} );

    const unsigned int k = source_range.max();

    assert( source_range.min() == 1 );
    assert( source_range.max() <= n );
    assert( sigma.isstrictlyascending() );

    std::vector<int> rho_values( n+1-k, -1 );
    
    int index = 0;
    const std::vector<int>& sigma_values = sigma.getvalues();

    assert( sigma_values.size() == k );
    
    for( auto value : IndexRange(0,n) )
    {
        int c = 0; 
        
        while( c < k and sigma_values[c] != value ) c++;
        
        if( c == k ) {
            rho_values[index] = value;
            index++;
        }
    }

    Assert( index == n-k+1, index, n, k, n-k+1 );

    IndexMap rho( IndexRange(0,n-k), IndexRange(0,n), rho_values );

    for( auto value : IndexRange(1,  k) ) assert( not   rho.has_value_in_range( sigma[value] ) );
    for( auto value : IndexRange(0,n-k) ) assert( not sigma.has_value_in_range(   rho[value] ) );

    return rho;
}




int sign_of_rho_sigma( const IndexMap& sigma )
{
    const auto source_range = sigma.getSourceRange();
    const auto target_range = sigma.getTargetRange();

    assert( not target_range.isempty() );
    
    const unsigned int n = target_range.max();

    assert( target_range.min() == 0 );
    assert( target_range.max() == n );

    if( source_range.isempty() )
        return 1;

    const unsigned int k = source_range.max();

    assert( source_range.min() == 1 );
    assert( source_range.max() <= n );
    assert( sigma.isstrictlyascending() );

    const auto rho = complement_sigma( sigma, n );

    assert( not rho.getSourceRange().isempty() and rho.getSourceRange().max() == n-k );

    const auto sigma_values = sigma.getvalues();
    const auto rho_values   =   rho.getvalues();

    auto values = std::vector<int>();
    values.reserve(n+1);

    for( auto value :   rho_values ) values.push_back( value );
    for( auto value : sigma_values ) values.push_back( value );
    
    int swaps = 0;
    
    for( int i =   0; i < values.size(); ++i )
    for( int j = i+1; j < values.size(); ++j )
    {
        if( i != j ) assert( values[i] != values[j] );

        if( values[i] > values[j] ) 
        {
            std::swap( values[i], values[j] );
            swaps++;
        }
    }

    for( int i = 1; i < values.size(); ++i )
        assert( values[i-1] < values[i] );
    
    return (swaps % 2 == 0) ? 1 : -1;
}

// IndexMap IndexMap::skip( int i ) const 
// {
//     check();
//     assert( src.contains(i) );
//     IndexMap ret( *this );
//     ret.src = IndexRange( src.min(), src.max() - 1 );
//     ret.values.erase( ret.values.begin() + ( i - src.min() ) );
//     ret.check();
//     return ret;
// }

// IndexMap IndexMap::shiftup() const 
// {
//     check();
//     IndexMap ret( *this );
//     ret.src = IndexRange( src.min()+1, src.max()+1 );
//     ret.check();
//     return ret;
// }

// IndexMap IndexMap::attachbefore( int to ) const 
// {
//     check();
//     assert( dest.contains(to) );
//     IndexMap ret( *this );
//     ret.src = IndexRange( src.min() - 1, src.max() );
//     ret.values.insert( ret.values.begin(), to );
//     ret.check();
//     return ret;
// }


