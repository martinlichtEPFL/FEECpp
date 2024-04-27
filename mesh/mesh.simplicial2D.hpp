#ifndef INCLUDEGUARD_MESH_SIMPLICIAL_2D_HPP
#define INCLUDEGUARD_MESH_SIMPLICIAL_2D_HPP


#include <array>
#include <string>
#include <utility>
#include <vector>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"
#include "coordinates.hpp"
#include "mesh.hpp"


/*******************
****  
****  
****  MeshSimplicial2D Class 
****  
****  - specialized mesh class for finite two-dimensional simplicial complexes
****    
****    
*******************/


class MeshSimplicial2D
: public Mesh
{

    public:
    
        /* Constructors */
        
        explicit MeshSimplicial2D( int outerdim = 2 );
        
        MeshSimplicial2D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,3>>& triangle_vertices
        );
        
        MeshSimplicial2D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,3>>& triangle_edges,
            const std::vector<int              >& edge_firstparent_triangle,
            const std::vector<std::array<int,3>>& triangle_nextparents_of_edges,
            const std::vector<std::array<int,3>>& triangle_vertices,
            const std::vector<int              >& vertex_firstparent_triangle,
            const std::vector<std::array<int,3>>& triangle_nextparents_of_vertices,
            const std::vector<std::array<int,2>>& edge_vertices,
            const std::vector<int              >& vertex_firstparent_edge,
            const std::vector<std::array<int,2>>& edge_nextparents_of_vertices
        );
        
        
        /* standard methods for operators */
        
        MeshSimplicial2D( const MeshSimplicial2D& ) = default;
        MeshSimplicial2D& operator=( const MeshSimplicial2D& ) = default;
        MeshSimplicial2D( MeshSimplicial2D&& ) = default;
        MeshSimplicial2D& operator=( MeshSimplicial2D&& ) = default;
        virtual ~MeshSimplicial2D();
        
        /* standard interface */
        
        virtual void check() const;
        
        // virtual void print( std::ostream& out ) const override;

        virtual std::string text() const override;

        /* OTHER METHODS */
        
        
        bool is_equal_to( const MeshSimplicial2D& mesh ) const;
        
        /* inherited methods */
        
        virtual bool has_dimension_counted( int dim ) const override;
        
        virtual int count_simplices( int dim ) const override;
        
        virtual bool has_subsimplices_listed( int sup, int sub ) const override;
        
        virtual IndexMap getsubsimplices( int sup, int sub, int cell ) const override;
        
        virtual bool has_supersimplices_listed( int sup, int sub ) const override;
        
        virtual const std::vector<int> getsupersimplices( int sup, int sub, int cell ) const override;
        
        
        virtual SimplexFlag get_flag( int dim, int index ) const override;
        
        virtual void set_flag( int dim, int index, SimplexFlag flag ) override;
        
        
        
        /* General management */
        
        
        /* count the simplices of a certain type */
        int count_triangles() const;
        
        int count_edges()     const;
        
        int count_vertices()  const;
        
        
        /* subsimplex relation of triangles and edges */
        
        bool contains_triangle_edge( int t, int e ) const;
        
        int indexof_triangle_edge( int t, int e ) const;
        
        int get_triangle_edge( int t, int ei ) const;
        
        const std::array<int,3> get_triangle_edges( int t ) const;
        
        
        /* subsimplex relation of triangles and vertices */
        
        bool contains_triangle_vertex( int t, int v ) const;
        
        int indexof_triangle_vertex( int t, int v ) const;
        
        int get_triangle_vertex( int t, int vi ) const;
        
        const std::array<int,3> get_triangle_vertices ( int t ) const;
        
        
        /* subsimplex relation of edges and vertices */
        
        bool contains_edge_vertex( int e, int v ) const;
        
        int indexof_edge_vertex( int e, int v ) const;
        
        int get_edge_vertex( int e, int vi ) const;
        
        const std::array<int,2> get_edge_vertices ( int e ) const;
        
        
        /* triangle parents of an edge */
        
        int count_edge_triangle_parents( int e ) const;
        
        int get_edge_firstparent_triangle( int e ) const;
        
        int get_edge_nextparent_triangle( int e, int t ) const;
        
        int get_triangle_nextparent_of_edge( int t, int ei ) const;
        
        bool is_triangle_edge_parent( int t, int e ) const;
        
        int indexof_triangle_edge_parent( int t, int e ) const;
        
        std::vector<int> get_triangle_parents_of_edge( int e ) const;
        
        
        /* triangle parents of a vertex */
        
        int count_vertex_triangle_parents( int v ) const;
        
        int get_vertex_firstparent_triangle( int v ) const;
        
        int get_vertex_nextparent_triangle( int v, int t ) const;
        
        int get_triangle_nextparent_of_vertex( int t, int vi ) const;
        
        bool is_triangle_vertex_parent( int t, int v ) const;
        
        int indexof_triangle_vertex_parent( int t, int v ) const;
        
        std::vector<int> get_triangle_parents_of_vertex( int v ) const;
        
        
        /* edge parents of a vertex */
        
        int count_vertex_edge_parents( int v ) const;
        
        int get_vertex_firstparent_edge( int v ) const;
        
        int get_vertex_nextparent_edge( int v, int e ) const;
        
        int get_edge_nextparent_of_vertex( int e, int vi ) const;
        
        bool is_edge_vertex_parent( int e, int v ) const;
        
        int indexof_edge_vertex_parent( int e, int v ) const;
        
        std::vector<int> get_edge_parents_of_vertex( int v ) const;
        
        
        /* refinement */
        
        void bisect_edge( int e );

        // void longest_edge_bisection( std::vector<int> edges );
        void longest_edge_bisection_recursive( const std::vector<int>& edges );
        void longest_edge_bisection_recursive( const int e );

        // void newest_vertex_bisection( std::vector<int> edges );
        void newest_vertex_bisection_recursive( const std::vector<int>& edges );
        void newest_vertex_bisection_recursive( int e );
        
        void uniformrefinement();
        
        void uniformrefinement( int levels );
        
        void midpoint_refinement( int t );
        
        void midpoint_refinement_global();
        
        
        /* other things */
        
        FloatVector get_triangle_midpoint( int t ) const;
        FloatVector get_edge_midpoint( int e ) const;
        Float get_edge_length( int e ) const;
        
        int get_oldest_edge( int t ) const;
                
        void merge( const MeshSimplicial2D& );
        
        
        /* TikZ */
        
        std::string outputTikZ() const;

        std::string outputSVG( 
            Float stroke_width = 0.01,
            const std::string& fill   = "white",
            const std::string& stroke = "black",
            const FloatVector* triangle_red   = nullptr,
            const FloatVector* triangle_green = nullptr,
            const FloatVector* triangle_blue  = nullptr
        ) const;

        std::string outputLinearSVG( 
            const FloatVector& triangle_red,
            const FloatVector& triangle_green,
            const FloatVector& triangle_blue, 
            Float stroke_width = 0.01,
            const std::string& fill   = "red",
            const std::string& stroke = "blue"
        ) const;


        /* other */ 
        
        virtual std::size_t memorysize() const override;
        
        
        
    private:

        int counter_triangles;
        int counter_edges;
        int counter_vertices;
        
        std::vector< std::array<int,3> > data_triangle_edges;
        std::vector< int               > data_edge_firstparent_triangle;
        std::vector< std::array<int,3> > data_triangle_nextparents_of_edges;
        
        std::vector< std::array<int,3> > data_triangle_vertices;
        std::vector< int               > data_vertex_firstparent_triangle;
        std::vector< std::array<int,3> > data_triangle_nextparents_of_vertices;
        
        std::vector< std::array<int,2> > data_edge_vertices;
        std::vector< int               > data_vertex_firstparent_edge;
        std::vector< std::array<int,2> > data_edge_nextparents_of_vertices;
        
        
        std::vector<SimplexFlag> flags_triangles;
        std::vector<SimplexFlag> flags_edges;
        std::vector<SimplexFlag> flags_vertices;


    public:

        inline bool operator==( const MeshSimplicial2D& m2 ) const 
        {
            return this->is_equal_to( m2 );
        }

        inline bool operator!=( const MeshSimplicial2D& m2 ) const 
        {
            return !( *this == m2 );
        }



};




// inline std::ostream& operator<<( std::ostream& os, const MeshSimplicial2D& mt2d )
// {
//     mt2d.print( os );
//     return os;
// }










#endif
