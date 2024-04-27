#ifndef INCLUDEGUARD_COMBINATORICS_MULTIINDEX_HPP
#define INCLUDEGUARD_COMBINATORICS_MULTIINDEX_HPP


#include <functional>
#include <initializer_list>
// #include <ostream>
#include <string>
#include <vector>

#include "../basic.hpp"

#include "indexrange.hpp"
#include "indexmap.hpp"


/********************
**** 
****  Class that describes multiindices over an IndexRange 
****  basic arithmetic functionality
**** 
********************/

class MultiIndex final
: public IndexMap
{

    public:

        /* Constructors */
        
        explicit MultiIndex( const IndexRange& ir );
        MultiIndex( const IndexRange& ir, const std::vector<int>& );
        MultiIndex( const IndexRange&, const std::function<int(int)>& );
        MultiIndex( const IndexRange&, const std::initializer_list<int>& );
        
        /* standard interface */ 
        
        MultiIndex( const MultiIndex& )             = default;
        MultiIndex& operator =( const MultiIndex& ) = default;
        MultiIndex( MultiIndex&& )                  = default;
        MultiIndex& operator =( MultiIndex&& )      = default;
        
        virtual ~MultiIndex() = default;
        
        /* standard methods */
        
        void check() const;

        std::string text( bool embellish = false ) const;
        
        // void print( std::ostream&, bool embellish = false ) const;

        // void lg() const { LOG << text() << nl; }
        
        
        
        /* OTHER METHODS */
        
        IndexRange getIndexRange() const;

        const std::vector<int>& getvalues() const;

                

        int absolute() const;

        int factorial() const;
        
        Float factorial_numerical() const;
        
        int min() const;

        int max() const;
        

        void add( int p );

        void sub( int p );
        
        void add( int p, int n );

        void sub( int p, int n );
        
        void add( const MultiIndex& mi );

        void sub( const MultiIndex& mi );
                

        bool is_comparable_with( const MultiIndex& mi ) const;

        bool is_equal_to( const MultiIndex& mi ) const;
                
        bool is_less_than( const MultiIndex& mi ) const;

};



inline MultiIndex ZeroMultiIndex( const IndexRange& ir )
{
    return MultiIndex( ir, std::vector<int>( ir.cardinality(), 0 ) );
}
        




template<typename Stream>
inline Stream& operator<<( Stream&& os, const MultiIndex& mi )
{
    mi.check();
    os << mi.text(); // mi.print( os );
    return os;
}

inline MultiIndex& operator+=( MultiIndex& left, int right )
{
    left.check();
    left.add( right );
    return left;
}

inline MultiIndex& operator-=( MultiIndex& left, int right )
{
    left.check();
    left.sub( right );
    return left;
}

inline MultiIndex& operator+=( MultiIndex& left, const MultiIndex& right )
{
    left.check();
    right.check();
    left.add( right );
    return left;
}

inline MultiIndex& operator-=( MultiIndex& left, const MultiIndex& right )
{
    left.check();
    right.check();
    left.sub( right );
    return left;
}


inline MultiIndex operator+( const MultiIndex& left, int right )
{
    left.check();
    MultiIndex ret = left;
    ret += right;
    ret.check();
    return ret;
}

inline MultiIndex operator-( const MultiIndex& left, int right )
{
    left.check();
    MultiIndex ret = left;
    ret -= right;
    ret.check();
    return ret;
}

inline MultiIndex operator+( const MultiIndex& left, const MultiIndex& right )
{
    left.check();
    right.check();
    MultiIndex ret = left;
    ret += right;
    ret.check();
    return ret;
}

inline MultiIndex operator-( const MultiIndex& left, const MultiIndex& right )
{
    left.check();
    right.check();
    MultiIndex ret = left;
    ret -= right;
    ret.check();
    return ret;
}


inline bool operator==( const MultiIndex& it, const MultiIndex& mi)
{
    it.check();
    mi.check();
    return it.is_equal_to( mi );
}
                
inline bool operator!=( const MultiIndex& it, const MultiIndex& mi)
{
    it.check();
    mi.check();
    return ! ( it == mi );
}

inline bool operator<( const MultiIndex& it, const MultiIndex& mi)
{
    it.check();
    mi.check();
    return it.is_less_than( mi );
}
                
inline bool operator>( const MultiIndex& it, const MultiIndex& mi)
{
    it.check();
    mi.check();
    return mi < it;
}
                
inline bool operator<=( const MultiIndex& it, const MultiIndex& mi)
{
    it.check();
    mi.check();
    return it < mi || it == mi;
}
                
inline bool operator>=( const MultiIndex& it, const MultiIndex& mi)
{
    it.check();
    mi.check();
    return it > mi || it == mi;
}
        
        
inline int absolute( const MultiIndex& it )
{
    return it.absolute();
}


inline int factorial( const MultiIndex& it )
{
    return it.factorial();
}




                
#endif
