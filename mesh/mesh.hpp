#ifndef INCLUDEGUARD_MESH_MESH_HPP
#define INCLUDEGUARD_MESH_MESH_HPP


// #include <ostream>
#include <utility>
#include <vector>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "coordinates.hpp"



/*******************
****
****    Datatype for the Simplex Flags 
****    
****    To be used to indicate properties such as Dirichlet boundary conditions
****    
****
*******************/


enum class SimplexFlag : unsigned char
{
    SimplexFlagNull      = 0x1B,
    SimplexFlagInvalid   = 0x77,
    SimplexFlagDirichlet = 0xF1 
};

// const SimplexFlag SimplexFlagNull    = 0x1B;
// const SimplexFlag SimplexFlagInvalid = 0x77;

// const SimplexFlag SimplexFlagDirichlet = 0xF1;



/*******************
****  
****  
****  Mesh Class 
****  
****  - contains Coordinate Class 
****  - intrinsic and exterior dimension
****  - provides relation about sub and supersimplices
****  - the topological data reflect a simplicial complex  
****    
****  
*******************/

class Mesh
{
    
    public:
        
        /* Constructors */
        
        Mesh( int inner, int outer );
        
        /* standard interface */
        
        Mesh( const Mesh& ) = default;
        Mesh& operator=( const Mesh& ) = default;
        Mesh( Mesh&& ) = default;
        Mesh& operator=( Mesh&& ) = default;
        virtual ~Mesh();
        
        
        /* standard methods for operators */
        
        void check() const;
        
        // void print( std::ostream& out ) const;
        
        virtual std::string text() const = 0;
        
        // // void lg() const { LOG << *this << nl; };
        
        
        /* OTHER METHODS */
        
        static const int nullindex; 
        
        static int is_not_nullindex( int i ){ return i != nullindex; }
        
        static int is_nullindex( int i ){ return i == nullindex; }
        
        
        
        /* Basic data */
        
        int getinnerdimension() const;
        
        int getouterdimension() const;
        
        Coordinates& getcoordinates();
        
        const Coordinates& getcoordinates() const;
        
        
        
        /* Static auxiliary functions */
        
        int index_from_pair( int sup, int sub ) const;
        
        void index_to_pair( int index, int& sup, int& sub ) const;
        
        int count_subsimplices( int sup, int sub ) const;
        
        
        
        /* Counting simplices */
        
        virtual bool has_dimension_counted( int dim ) const = 0;
        
        virtual int count_simplices( int dim ) const = 0;
        
        
        
        /* 
         * Accessing subsimplices 
         *
         * - check whether subsimplices are listed at all 
         * - get the subsimplex list of a cell
         * - test for subsimplex relation 
         * - get local index of subsimplex 
         * - get listed subsimplex 
         * 
         */
        
        virtual bool has_subsimplices_listed( int sup, int sub ) const = 0;
        
        virtual IndexMap getsubsimplices( int sup, int sub, int cell ) const = 0;
        
        virtual bool is_subsimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_subsimplex_index( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_subsimplex( int sup, int sub, int cellsup, int localindex ) const;
        
        int get_opposite_subsimplex_index( int sup, int sub, int cellsup, int localindex ) const;
        
        
        
        /* 
         * Accessing supersimplices 
         * 
         * - check whether supersimplices are listed at all 
         * - get the supersimplex list of a cell
         * - test for supersimplex relation 
         * - get local index of supersimplex 
         * 
         */
        
        virtual bool has_supersimplices_listed( int sup, int sub ) const = 0;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const = 0;
        
        virtual bool is_supersimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_firstparent_of_subsimplex( int sup, int sub, int cellsub ) const;
        
        virtual int get_nextparent_of_subsimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_nextparent_by_localindex( int sup, int sub, int cellsup, int localindex ) const;
        
//         virtual IndexMap getnextparents( int sup, int sub, int cell ) const;
        
        virtual int get_index_of_supersimplex( int sup, int sub, int cellsup, int cellsub ) const;
        
        virtual int get_supersimplex_by_index( int sup, int sub, int cellsub, int parentindex ) const;

        
        /* 
         * 
         * Setting and getting flags 
         * 
         */
        
        virtual SimplexFlag get_flag( int dim, int index ) const = 0;
        
        virtual void set_flag( int dim, int index, SimplexFlag flag ) = 0;
        
        void set_flags( int dim, SimplexFlag flag );
        
        const std::vector<SimplexFlag> get_flags( int dim ) const;
        
        void set_flags( int dim, std::vector<SimplexFlag> flags );
        
        void automatic_dirichlet_flags();
        
        void check_dirichlet_flags( bool check_for_full_dirichlet = true );
        
        
        /* 
         * Accessing geometric information about the mesh
         * 
         */
        
        Float getDiameter( int dim, int index ) const;
        
        Float getMaximumDiameter() const;
        Float getMinimumDiameter() const;
        
        Float getMeasure( int dim, int index ) const;
        
        Float getHeight( int dim, int index, int vertexindex ) const;
        
        Float getHeightQuotient( int dim, int index ) const;
        Float getHeightQuotient( int dim ) const;
        Float getHeightQuotient() const;
        
        Float getShapemeasure( int dim, int index ) const;
        Float getShapemeasure( int dim ) const;
        Float getShapemeasure() const;

        int getVertexPatchSize() const;
        
        int getSupersimplexSize( int dim ) const;
        
        Float getComparisonQuotient() const;

        Float getRadiiQuotient( int dim ) const;
        
        FloatVector get_midpoint( int dim, int index ) const;

        FloatVector get_random_point( int dim, int index ) const;

        FloatVector getPointFromBarycentric( int dim, int index, const FloatVector& barycoords ) const;

        
        int get_longest_edge_index( int dim, int index ) const;
        
        DenseMatrix getVertexCoordinateMatrix( int dim, int index ) const;
        
        DenseMatrix getTransformationJacobian( int dim, int index ) const;

        DenseMatrix getBarycentricProjectionMatrix( int dim, int index ) const; 
        
        DenseMatrix getGradientMatrix( int dim, int index ) const;
        
        DenseMatrix getGradientProductMatrix( int dim, int index ) const;
        
        DenseMatrix getGradientProductMatrixRightFactor( int dim, int index ) const;


        FloatVector transform_whitney_to_euclidean( int dim, const FloatVector& whitneyvalues, int zero_padding = 0 ) const;

        
        /* 
         * Manipulation
         * 
         */
        
        void shake_interior_vertices( Float intensity = 0.10, Float probability = 0.5 );
        

        virtual std::size_t memorysize() const = 0;
        
    private:
        
        int innerdimension;
        int outerdimension;
        
        Coordinates coordinates;
        
    private:
        
        // std::vector< std::vector< std::vector<IndexMap> > > auxdata;
        
    public: 
        
#if __cplusplus >= 201402L
        // const auto& getauxdata() { return auxdata; }
#else
        // const decltype(auxdata)& getauxdata() { return auxdata; }
#endif
        
        template<typename Stream>
        friend inline Stream& operator<<( Stream&& os, const Mesh& mesh )
        {
            os << mesh.text(); // mesh.print( os );
            return os;
        }

};







FloatVector get_random_barycentric_coordinates( int dim );

// static inline int countsubsimplices( int n, int k )
// {
//     assert( 0 <= k && k <= n );
//     return binomial_integer( n+1, k+1 );
// }











#endif
