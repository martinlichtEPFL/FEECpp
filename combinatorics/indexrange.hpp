#ifndef INCLUDEGUARD_COMBINATORICS_INDEXRANGE_HPP
#define INCLUDEGUARD_COMBINATORICS_INDEXRANGE_HPP


// #include <ostream>
#include <limits>
#include <string>

#include "../basic.hpp"

/*****
**
**  This class models a range of indices, i.e. integers,
**  of the form \{ min, min+1, ..., max \}.
**  The index range may be empty.
**
******/


class IndexRange final
{
    
    public:
        
        /* Constructors */
        
        IndexRange( int from, int to );
        
        IndexRange( const IndexRange& )             = default;
        IndexRange& operator =( const IndexRange& ) = default;
        IndexRange( IndexRange&& )                  = default;
        IndexRange& operator =( IndexRange&& )      = default;
        ~IndexRange()                               = default;
        
        /* standard methods */

        void check() const;
        
        std::string text( bool embellish = false ) const;
        
        // void print( std::ostream&, bool embellish = true ) const;

        // void lg() const { LOG << text() << nl; };
        
        /* OTHER METHODS */

        int min() const;
        
        int max() const;
        
        int cardinality() const;
        
        bool isempty() const;
        
        bool contains( int element ) const;
        
        bool contains( const IndexRange& subir ) const;
        
        bool isequal( const IndexRange& ir ) const;
        
        int element2position( int element ) const;
        
        int position2element( int position ) const;
        
        /* For each semantics */ 
        
        class ConstIterator {
            
            friend IndexRange;
        
            private: 
            
                int value;
                int minimum;
                int maximum;
                bool is_end;
                
                explicit ConstIterator( int value, int minimum, int maximum, bool is_end ) : value(value), minimum(minimum), maximum(maximum), is_end(is_end)
                { }
                
            public:
                
                inline int operator*() const
                {
                    assert( !is_end );
                    if( minimum <= maximum ) assert( minimum <= value && value <= maximum );
                    return value;                
                }
                
                inline ConstIterator operator++()
                {
                    assert( !is_end );
                    if( value == maximum )
                        is_end = true;
                    else 
                        ++value; 
                    return *this;
                } // pre-increment
                
                inline ConstIterator operator++( int )
                {
                    assert( !is_end );
                    auto ret = *this;
                    ++(*this); 
                    return ret;
                } // post-increment
                
                inline bool operator!=( const ConstIterator& irit ) const 
                { 
                    return ( is_end != irit.is_end ) || ( !is_end && value != irit.value );
                }
                    
                inline bool operator==( const ConstIterator& irit ) const 
                { 
                    return !( *this != irit );
                }
                    
        };
        
        inline ConstIterator begin() const { return ConstIterator(minimum,minimum,maximum,minimum>maximum); }
        
        inline ConstIterator end()   const { return ConstIterator(minimum,minimum,maximum,true ); }
        
        
    private:

        int minimum;
        int maximum;
        
};

template<typename Stream>
inline Stream& operator<<( Stream&& os, const IndexRange& ir )
{
    os << ir.text(); // ir.print( os );
    return os;
}

inline bool operator== ( const IndexRange& ir1, const IndexRange& ir2 )
{
    return ir1.isequal( ir2 );
}

inline bool operator!= ( const IndexRange& ir1, const IndexRange& ir2 )
{
    return !( ir1 == ir2 );
}


static const IndexRange  NonNegativeIntegers = IndexRange( 0, std::numeric_limits<int>::max() );

static const IndexRange  PositiveIntegers    = IndexRange( 1, std::numeric_limits<int>::max() );



#endif
