#ifndef INCLUDEGUARD_MESH_SIMPLICIAL_1D_HPP
#define INCLUDEGUARD_MESH_SIMPLICIAL_1D_HPP


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
****  MeshSimplicial1D Class 
****  
****  - specialized mesh class for finite one-dimensional simplicial complexes 
****  - looks like an undirected graph without loops or parallel edges.
****  
****    Content:
****    - for every edge, we save the 2 vertices 
****    - for every vertex, we save the first parent edge
****    - for every edge, we save the next parents of each vertex. 
****    
****    
*******************/


class MeshSimplicial1D
: public Mesh
{

    public:
    
        /* Constructors */
        
        explicit MeshSimplicial1D( int outerdim = 1 );
        
        MeshSimplicial1D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,2>>& edge_vertices
        );
        
        MeshSimplicial1D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,2>>& edge_vertices,
            const std::vector<std::array<int,2>>& edge_nextparents_of_vertices,
            const std::vector<int              >& vertex_firstparent_edge
        );
        
        /* standard methods for operators */
        
        MeshSimplicial1D( const MeshSimplicial1D& ) = default;
        MeshSimplicial1D& operator=( const MeshSimplicial1D& ) = default;
        MeshSimplicial1D( MeshSimplicial1D&& ) = default;
        MeshSimplicial1D& operator=( MeshSimplicial1D&& ) = default;
        virtual ~MeshSimplicial1D();
        
        /* standard interface */
        
        virtual void check() const;
        
        // virtual void print( std::ostream& out ) const override;
        
        virtual std::string text() const override;
        
        
        /* OTHER METHODS */
        
        bool is_equal_to( const MeshSimplicial1D& mesh ) const;
        
        
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
        int count_edges()    const;
        int count_vertices() const;
        
        
        /* subsimplex relation of edges and vertices */
        
        bool contains_edge_vertex( int e, int v ) const;
        
        int indexof_edge_vertex( int e, int v ) const;
        
        int get_edge_vertex( int e, int vi ) const;
        
        const std::array<int,2> get_edge_vertices ( int e ) const;
        
        
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
        
        void uniformrefinement();
        
        void improved_uniformrefinement();
        
        void midpoint_refinement( int e );
        
        void midpoint_refinement_global();
        
        
        
        /* other things */
        
        FloatVector get_edge_midpoint( int e ) const;
        
        void merge( const MeshSimplicial1D& other );

        virtual std::size_t memorysize() const override;
        
    private:

        int counter_edges;
        int counter_vertices;
        
        std::vector< std::array<int,2> > data_edge_vertices;
        std::vector< int               > data_vertex_firstparent_edge;
        std::vector< std::array<int,2> > data_edge_nextparents_of_vertices;
        
        std::vector<SimplexFlag> flags_edges;
        std::vector<SimplexFlag> flags_vertices;


    public:

        inline bool operator==( const MeshSimplicial1D& m2 ) const 
        {
            return this->is_equal_to( m2 );
        }

        inline bool operator!=( const MeshSimplicial1D& m2 ) const 
        {
            return !( *this == m2 );
        }

        
        
};




// inline std::ostream& operator<<( std::ostream& os, const MeshSimplicial1D& mt1d )
// {
//     mt1d.print( os );
//     return os;
// }







#endif
