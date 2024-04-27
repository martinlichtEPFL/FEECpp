#ifndef INCLUDEGUARD_MESH_COORDINATES_HPP
#define INCLUDEGUARD_MESH_COORDINATES_HPP


#include <array>
// #include <ostream>
#include <vector>

#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../dense/densematrix.hpp"



/*******************
****  
****  Class for Coordiante Collections 
****  
****  - Most basic functionality
****  
****  
****  
****  
*******************/




class Coordinates
{

    public:

        Coordinates( const Coordinates& ) = default;
        Coordinates& operator=( const Coordinates& ) = default;
        Coordinates( Coordinates&& ) = default;
        Coordinates& operator=( Coordinates&& ) = default;

        
        
        Coordinates( int dimension, int number );
        Coordinates( int dimension, int number, const std::vector<Float>& );
        virtual ~Coordinates();
        
        void check() const;
        // void print( std::ostream& ) const;
        std::string text() const;
        
        // // void lg() const { LOG << *this << nl; };
        

        // void read( std::istream& ) ;

        int getdimension() const;
        int getnumber() const;
        IndexRange getIndexRange() const;
        
        /* get/set coordinates as per point  */
        
        Float getdata( int n, int d ) const;
        void setdata( int n, int d, Float v );
        
        /* get range of coordinates */
        
        Float getmin( int d) const;
        Float getmax( int d) const;
        
        /* get/set points as vectors  */
        
        FloatVector getvectorclone( int n ) const;
        FloatVector getvectorclone( int n, Float s ) const;
        void loadvector( int n, const FloatVector& value );
        void loadvector( int n, const FloatVector& value, Float s );
        
        /* get/set coordinates as vectors  */
        
        FloatVector getdimensionclone( int d, Float s = 1.0 ) const;
        void loaddimension( int d, const FloatVector& value, Float s = 1.0 );
        
        /* transform all coordinates  */
        
        void scale( Float alpha );
        void scale( FloatVector alphas );
        void shift( const FloatVector& add );
        void lineartransform( const LinearOperator& op );
        
        /* Add additional coordiantes */
        
        void append( const Coordinates& co );
        void append( const FloatVector& v );
        // void append( const std::vector<FloatVector>& );
        
        void addcapacity( int additional_capacity );
        void addcoordinates( int add_number );
        
        /* Obtain information about reference transformation of simplex */
        
        DenseMatrix getLinearPart( const IndexMap& im ) const;
        FloatVector getShiftPart( const IndexMap& im ) const;
        
        FloatVector getCenter() const;


        /* other */ 

        std::vector<Float>& raw();
        const std::vector<Float>& raw() const;
        
        std::size_t memorysize() const;
        
    private:
            
        int dimension;
        int number;
        std::vector<Float> data;

    public:

        static bool is_equal_to( const Coordinates& coords_left, const Coordinates& coords_right );

        friend inline bool operator==( const Coordinates& coords_left, const Coordinates& coords_right )
        {
            return is_equal_to( coords_left, coords_right );
        }

        friend inline bool operator!=( const Coordinates& coords_left, const Coordinates& coords_right )
        {
            return ! ( coords_left == coords_right );
        }

        friend inline std::ostream& operator<<( std::ostream& os, const Coordinates& co )
        {
            os << co.text(); // co.print( os );
            return os;
        }


};







#endif
