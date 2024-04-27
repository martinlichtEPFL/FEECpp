#ifndef INCLUDEGUARD_COMBINATORICS_INDEXMAP_HPP
#define INCLUDEGUARD_COMBINATORICS_INDEXMAP_HPP


#include <functional>
#include <initializer_list>
// #include <ostream>
#include <string>
#include <vector>

#include "../basic.hpp"

#include "indexrange.hpp"


/*****
**
**  This class models a mapping of indices.
**  It has a 'source' and 'target' range.
**  The mapping may be empty.
**
******/



class IndexMap
{
    
    public:
        
        /* Constructors */
        
        IndexMap( const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const std::initializer_list<int>& );
        
        IndexMap( const IndexRange&, const IndexRange&, const int& );
        IndexMap( const IndexRange&, const IndexRange&, const std::vector<int>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::function<int(int)>& );
        IndexMap( const IndexRange&, const IndexRange&, const std::initializer_list<int>& );
        
        /* standard interface */ 
        
        IndexMap()                              = delete;
        IndexMap( const IndexMap& )             = default;
        IndexMap& operator =( const IndexMap& ) = default;
        IndexMap( IndexMap&& )                  = default;
        IndexMap& operator =( IndexMap&& )      = default;
        virtual ~IndexMap()                     = default; // dtor virtualization is a bit shady here
        
        /* standard methods */
        
        void check() const;
        
        std::string text( bool embellish = true ) const;
        
        // void print( std::ostream&, bool embellish = true ) const;

        // void lg() const { LOG << text() << nl; };

        
        /* OTHER METHODS */
        
        const IndexRange& getSourceRange() const;
        
        const IndexRange& getTargetRange() const;
        
        bool isempty() const;
        
        bool isinjective() const;
        
        bool issurjective() const;
        
        bool isbijective() const;
        
        bool isstrictlyascending() const;
        
        
        
        
        int& at( int i ) &;
        
        const int& at( int i ) const &;
        
        int& operator[]( int i ) &;
        
        const int& operator[]( int i ) const &;
        
        const std::vector<int>& getvalues() const &;
                
        bool has_value_in_range( int value ) const;
        
        int preimageof( int value ) const;
        
        
        
        // IndexMap skip( int i ) const;
        // IndexMap shiftup() const;



        bool is_comparable_with( const IndexMap& im ) const;
        
        bool is_equal_to( const IndexMap& im ) const;
        
        bool is_less_than( const IndexMap& im ) const;
    
    private:
        
        IndexRange src;
        IndexRange dest;
        
        std::vector<int> values;
        
};



IndexMap mergeSigmas( const IndexMap& left, const IndexMap& right, int& sign );



inline IndexMap operator*( const IndexMap& leave, const IndexMap& enter )
{
    leave.check();
    enter.check();
    assert( enter.getTargetRange() == leave.getSourceRange() );

    IndexRange src  = enter.getSourceRange();
    IndexRange dest = leave.getTargetRange();
    IndexMap ret( src, dest, [ &leave, &enter ]( int i ) -> int { return leave[ enter[i] ]; } );

    ret.check();
    return ret;
}

inline bool operator==( const IndexMap& left, const IndexMap& right )
{
    left.check();
    right.check();
    assert( left.is_comparable_with( right ) );

    return left.is_equal_to( right );
}

inline bool operator!=( const IndexMap& left, const IndexMap& right )
{
    left.check();
    right.check();
    assert( left.is_comparable_with( right ) );

    return !( left.is_equal_to( right ) );
}

inline bool operator<( const IndexMap& left, const IndexMap& right )
{
    left.check();
    right.check();
    assert( left.is_comparable_with( right ) );

    return left.is_less_than( right );
}

template<typename Stream>
inline Stream& operator<<( Stream&& os, const IndexMap& im )
{
    im.check();
    os << im.text(); // im.print( os );
    return os;
}



inline IndexMap identityIndexMap( const IndexRange& ir )
{
    ir.check();

    IndexMap im( ir, ir, []( int i ) -> int { return i; } );
    im.check();

    return im;
}

inline IndexMap identityIndexMap( int low, int high )
{
    return identityIndexMap( IndexRange( low, high ) );
}



// inline IndexMap toss_entry( const IndexMap& original, int p )
// {
    
//     assert( not original.getSourceRange().isempty() );
//     assert( original.getSourceRange().contains(p) );
    
//     IndexRange new_source_range = IndexRange( original.getSourceRange().min()+1, original.getSourceRange().max() );

//     IndexMap im( new_source_range, original.getTargetRange(), [p,original]( int i ) -> int { 
//         if( i <= p ) 
//             return original.at(p-1);
//         else 
//             return original.at(p);
//     } );

//     // IndexMap im( new_source_range, original.getTargetRange(), 0 );
//     // for( int j = 1; j <= p; j++ ) im[j] = original[j-1];
//     // for( int j = p+1; j < original.getSourceRange().max(); j++ ) im[j] = original[j];
    
//     return im;
// }





IndexMap expand_zero( const IndexMap& im, int p );

IndexMap expand_one( const IndexMap& im, int p );

IndexMap complement_sigma( const IndexMap& sigma, const int n );

int sign_of_rho_sigma( const IndexMap& sigma );



// inline int fehlstelle( const IndexMap& sub, const IndexMap& super )
// {
//     sub.check();
//     super.check();
//     assert( sub.getSourceRange().getlength() == super.getSourceRange().getlength() + 1 );
//     for( int i : sub.getSourceRange() )
//         assert( super.has_value_in_range(i) );
//     for( int j : super.getSourceRange() )
//         if( ! sub.has_value_in_range(j) )
//             return j;
//         else 
//             unreachable();
// }



#endif
