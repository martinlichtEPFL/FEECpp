#ifndef INCLUDEGUARD_MESH_SIMPLICIAL_3D_HPP
#define INCLUDEGUARD_MESH_SIMPLICIAL_3D_HPP


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
****  MeshSimplicial3D Class 
****  
****  - specialized mesh class for finite two-dimensional simplicial complexes
****    
****    
*******************/


class MeshSimplicial3D
: public Mesh
{

    public:
    
        /* Constructors */
        
        explicit MeshSimplicial3D( int outerdim = 3 );
        
        MeshSimplicial3D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,4>>& tetrahedron_vertices
        );
        
        MeshSimplicial3D( 
            int outerdim,
            const Coordinates& coords,
            const std::vector<std::array<int,4>>& tetrahedron_faces,
            const std::vector<int              >& face_firstparent_tetrahedron,
            const std::vector<std::array<int,4>>& tetrahedron_nextparents_of_faces,
            const std::vector<std::array<int,6>>& tetrahedron_edges,
            const std::vector<int              >& edge_firstparent_tetrahedron,
            const std::vector<std::array<int,6>>& tetrahedron_nextparents_of_edges,
            const std::vector<std::array<int,4>>& tetrahedron_vertices,
            const std::vector<int              >& vertex_firstparent_tetrahedron,
            const std::vector<std::array<int,4>>& tetrahedron_nextparents_of_vertices,
            const std::vector<std::array<int,3>>& face_edges,
            const std::vector<int              >& edge_firstparent_face,
            const std::vector<std::array<int,3>>& face_nextparents_of_edges,
            const std::vector<std::array<int,3>>& face_vertices,
            const std::vector<int              >& vertex_firstparent_face,
            const std::vector<std::array<int,3>>& face_nextparents_of_vertices,
            const std::vector<std::array<int,2>>& edge_vertices,
            const std::vector<int              >& vertex_firstparent_edge,
            const std::vector<std::array<int,2>>& edge_nextparents_of_vertices
        );

        /* standard interface */
        
        MeshSimplicial3D( const MeshSimplicial3D& ) = default;
        MeshSimplicial3D& operator=( const MeshSimplicial3D& ) = default;
        MeshSimplicial3D( MeshSimplicial3D&& ) = default;
        MeshSimplicial3D& operator=( MeshSimplicial3D&& ) = default;
        virtual ~MeshSimplicial3D();
        
        /* standard methods for operators */
        
        virtual void check() const;
        
        // virtual void print( std::ostream& out ) const override;

        virtual std::string text() const override;


        /* OTHER METHODS */
        
        bool is_equal_to( const MeshSimplicial3D& mesh ) const;
        
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
        int count_tetrahedra() const;
        int count_faces()      const;
        int count_edges()      const;
        int count_vertices()   const;
        
        
        /* subsimplex relation of tetrahedra and faces */
        
        bool contains_tetrahedron_face( int t, int f ) const;
        
        int indexof_tetrahedron_face( int t, int f ) const;
        
        int get_tetrahedron_face( int t, int fi ) const;
        
        const std::array<int,4> get_tetrahedron_faces( int t ) const;
        
        
        /* subsimplex relation of tetrahedra and edges */
        
        bool contains_tetrahedron_edge( int t, int e ) const;
        
        int indexof_tetrahedron_edge( int t, int e ) const;
        
        int get_tetrahedron_edge( int t, int ei ) const;
        
        const std::array<int,6> get_tetrahedron_edges( int t ) const;
        
        
        /* subsimplex relation of tetrahedra and vertices */
        
        bool contains_tetrahedron_vertex( int t, int v ) const;
        
        int indexof_tetrahedron_vertex( int t, int v ) const;
        
        int get_tetrahedron_vertex( int t, int vi ) const;
        
        const std::array<int,4> get_tetrahedron_vertices ( int t ) const;
        
        
        
        /* subsimplex relation of faces and edges */
        
        bool contains_face_edge( int f, int e ) const;
        
        int indexof_face_edge( int f, int e ) const;
        
        int get_face_edge( int f, int ei ) const;
        
        const std::array<int,3> get_face_edges( int f ) const;
        
        
        /* subsimplex relation of faces and vertices */
        
        bool contains_face_vertex( int f, int v ) const;
        
        int indexof_face_vertex( int f, int v ) const;
        
        int get_face_vertex( int f, int vi ) const;
        
        const std::array<int,3> get_face_vertices ( int f ) const;
        
        
        /* subsimplex relation of edges and vertices */
        
        bool contains_edge_vertex( int e, int v ) const;
        
        int indexof_edge_vertex( int e, int v ) const;
        
        int get_edge_vertex( int e, int vi ) const;
        
        const std::array<int,2> get_edge_vertices ( int e ) const;
        
        
        
        
        
        
        
        
        
        /* tetrahedron parents of a face */
        
        int count_face_tetrahedron_parents( int f ) const;
        
        int get_face_firstparent_tetrahedron( int f ) const;
        
        int get_face_nextparent_tetrahedron( int f, int t ) const;
        
        int get_tetrahedron_nextparent_of_face( int t, int fi ) const;
        
        bool is_tetrahedron_face_parent( int t, int f ) const;
        
        int indexof_tetrahedron_face_parent( int t, int f ) const;
        
        std::vector<int> get_tetrahedron_parents_of_face( int f ) const;
        
        
        /* tetrahedron parents of an edge */
        
        int count_edge_tetrahedron_parents( int e ) const;
        
        int get_edge_firstparent_tetrahedron( int e ) const;
        
        int get_edge_nextparent_tetrahedron( int e, int t ) const;
        
        int get_tetrahedron_nextparent_of_edge( int t, int ei ) const;
        
        bool is_tetrahedron_edge_parent( int t, int e ) const;
        
        int indexof_tetrahedron_edge_parent( int t, int e ) const;
        
        std::vector<int> get_tetrahedron_parents_of_edge( int e ) const;
        
        
        /* tetrahedron parents of a vertex */
        
        int count_vertex_tetrahedron_parents( int v ) const;
        
        int get_vertex_firstparent_tetrahedron( int v ) const;
        
        int get_vertex_nextparent_tetrahedron( int v, int t ) const;
        
        int get_tetrahedron_nextparent_of_vertex( int t, int vi ) const;
        
        bool is_tetrahedron_vertex_parent( int t, int v ) const;
        
        int indexof_tetrahedron_vertex_parent( int t, int v ) const;
        
        std::vector<int> get_tetrahedron_parents_of_vertex( int v ) const;
        
        
        /* face parents of an edge */
        
        int count_edge_face_parents( int e ) const;
        
        int get_edge_firstparent_face( int e ) const;
        
        int get_edge_nextparent_face( int e, int f ) const;
        
        int get_face_nextparent_of_edge( int f, int ei ) const;
        
        bool is_face_edge_parent( int f, int e ) const;
        
        int indexof_face_edge_parent( int f, int e ) const;
        
        std::vector<int> get_face_parents_of_edge( int e ) const;
        
        
        /* face parents of a vertex */
        
        int count_vertex_face_parents( int v ) const;
        
        int get_vertex_firstparent_face( int v ) const;
        
        int get_vertex_nextparent_face( int v, int f ) const;
        
        int get_face_nextparent_of_vertex( int f, int vi ) const;
        
        bool is_face_vertex_parent( int f, int v ) const;
        
        int indexof_face_vertex_parent( int f, int v ) const;
        
        std::vector<int> get_face_parents_of_vertex( int v ) const;
        
        
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
        void longest_edge_bisection_recursive( int e );
        
        void uniformrefinement();
        
        void uniformrefinement( int levels );
        
        void midpoint_refinement( int t );
        
        void midpoint_refinement_global();
        
        
        /* other things */
        
        FloatVector get_tetrahedron_midpoint( int t ) const;
        FloatVector get_face_midpoint( int f ) const;
        FloatVector get_edge_midpoint( int e ) const;
        Float get_edge_length( int e ) const;
        
        void merge( const MeshSimplicial3D& other );
        

        int get_oldest_edge( int t ) const;
        


        /* TikZ */
        
        std::string outputTikZ( bool boundary_only = false ) const;


        /* other */ 
        
        virtual std::size_t memorysize() const override;
        
    private:

        int counter_tetrahedra;
        int counter_faces;
        int counter_edges;
        int counter_vertices;
        
        std::vector< std::array<int,4> > data_tetrahedron_faces;
        std::vector< int               > data_face_firstparent_tetrahedron;
        std::vector< std::array<int,4> > data_tetrahedron_nextparents_of_faces;
        
        std::vector< std::array<int,6> > data_tetrahedron_edges;
        std::vector< int               > data_edge_firstparent_tetrahedron;
        std::vector< std::array<int,6> > data_tetrahedron_nextparents_of_edges;
        
        std::vector< std::array<int,4> > data_tetrahedron_vertices;
        std::vector< int               > data_vertex_firstparent_tetrahedron;
        std::vector< std::array<int,4> > data_tetrahedron_nextparents_of_vertices;
        
        std::vector< std::array<int,3> > data_face_edges;
        std::vector< int               > data_edge_firstparent_face;
        std::vector< std::array<int,3> > data_face_nextparents_of_edges;
        
        std::vector< std::array<int,3> > data_face_vertices;
        std::vector< int               > data_vertex_firstparent_face;
        std::vector< std::array<int,3> > data_face_nextparents_of_vertices;
        
        std::vector< std::array<int,2> > data_edge_vertices;
        std::vector< int               > data_vertex_firstparent_edge;
        std::vector< std::array<int,2> > data_edge_nextparents_of_vertices;


        std::vector<SimplexFlag> flags_tetrahedra;
        std::vector<SimplexFlag> flags_faces;
        std::vector<SimplexFlag> flags_edges;
        std::vector<SimplexFlag> flags_vertices;


    public:

        inline bool operator==( const MeshSimplicial3D& m2 ) const 
        {
            return this->is_equal_to( m2 );
        }

        inline bool operator!=( const MeshSimplicial3D& m2 ) const 
        {
            return !( *this == m2 );
        }



};




// inline std::ostream& operator<<( std::ostream& os, const MeshSimplicial3D& mt3d )
// {
//     mt3d.print( os );
//     return os;
// }















#endif
