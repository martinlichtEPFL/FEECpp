
#include <cmath>
#include <algorithm>
#include <sstream>
#include <list>
#include <map>
#include <stack> 
#include <string>
#include <utility>
#include <vector>


#include "../basic.hpp"
#include "../utility/sorthack.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "mesh.hpp"

#include "mesh.simplicial2D.hpp"




MeshSimplicial2D::MeshSimplicial2D( int outerdim )
:
    Mesh( 2, outerdim ),
    
    counter_triangles(0),
    counter_edges(0),
    counter_vertices(0),
    
    data_triangle_edges(0),
    data_edge_firstparent_triangle(0),
    data_triangle_nextparents_of_edges(0),
    
    data_triangle_vertices(0),
    data_vertex_firstparent_triangle(0),
    data_triangle_nextparents_of_vertices(0),
    
    data_edge_vertices(0),
    data_vertex_firstparent_edge(0),
    data_edge_nextparents_of_vertices(0),
    
    flags_triangles( 0, SimplexFlag::SimplexFlagNull ),
    flags_edges    ( 0, SimplexFlag::SimplexFlagNull ),
    flags_vertices ( 0, SimplexFlag::SimplexFlagNull )
    
{
    MeshSimplicial2D::check();
}



MeshSimplicial2D::MeshSimplicial2D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,3>>& triangle_vertices
)
:
    Mesh( 2, outerdim ),
    
    counter_triangles( triangle_vertices.size() ),
    counter_edges( 0 ),
    counter_vertices( 0 ),
    
    data_triangle_edges( triangle_vertices.size(), { nullindex, nullindex, nullindex } ),
    data_edge_firstparent_triangle( 0 ),
    data_triangle_nextparents_of_edges( triangle_vertices.size(), { nullindex, nullindex, nullindex } ),
    
    data_triangle_vertices( triangle_vertices ),
    data_vertex_firstparent_triangle( 0 ),
    data_triangle_nextparents_of_vertices( triangle_vertices.size(), { nullindex, nullindex, nullindex } ),
    
    data_edge_vertices( 0 ),
    data_vertex_firstparent_edge( 0 ),
    data_edge_nextparents_of_vertices( 0 ),
    
    flags_triangles( counter_triangles, SimplexFlag::SimplexFlagNull ),
    flags_edges    ( 0, SimplexFlag::SimplexFlagNull ),
    flags_vertices ( 0, SimplexFlag::SimplexFlagNull )
    

{
    
    getcoordinates() = coords;
    
    /* 0. sort the triangle vertices */

    for( auto& tri : data_triangle_vertices )
      std::sort( tri.begin(), tri.end() );

    /* 0. sort the triangles altogether */
    {    
//       std::sort( data_triangle_vertices.begin(), data_triangle_vertices.end() ); 
      sorthack( data_triangle_vertices );
      auto it = std::unique( data_triangle_vertices.begin(), data_triangle_vertices.end() );
      data_triangle_vertices.resize( it - data_triangle_vertices.begin() );
    }

    
    /* 1. create all edges, allocate memory */
    
    data_edge_vertices.resize( counter_triangles * 3 );
    for( int t  = 0; t  < counter_triangles; t++  )
    for( int ei = 0; ei <                 3; ei++ )
    {
      data_edge_vertices[ 0 * counter_triangles + t ] = { data_triangle_vertices[t][0], data_triangle_vertices[t][1] };
      data_edge_vertices[ 1 * counter_triangles + t ] = { data_triangle_vertices[t][0], data_triangle_vertices[t][2] };
      data_edge_vertices[ 2 * counter_triangles + t ] = { data_triangle_vertices[t][1], data_triangle_vertices[t][2] };
    }
    
    
    /* Due to an obscure bug in the C++ standard itself (!), we cannot use the standard sort function here */
    /* For that reason, a hand-written sort is used. No complexity estimates are guaranteed. */
    
//     std::sort( data_edge_vertices.begin(), data_edge_vertices.end() ); 
    sorthack( data_edge_vertices );
    auto it = std::unique( data_edge_vertices.begin(), data_edge_vertices.end() );
    data_edge_vertices.resize( it - data_edge_vertices.begin() );
    
    counter_edges = data_edge_vertices.size();
    
    data_edge_firstparent_triangle.resize   ( counter_edges, nullindex );
    data_edge_nextparents_of_vertices.resize( counter_edges, { nullindex, nullindex } );
    
    /* 2. Count vertices, allocate memory */
    
    counter_vertices = 0;
    for( const auto& duple : data_edge_vertices )
    for( const int& vertex : duple )
      counter_vertices = counter_vertices < vertex ? vertex : counter_vertices; 
    counter_vertices += 1;
    
    data_vertex_firstparent_triangle.resize( counter_vertices, nullindex );
    data_vertex_firstparent_edge.resize( counter_vertices, nullindex );
    
      
    /* 3. For each vertex, set the first parent triangle and the neighboring parent triangles */
    
    for( int t =  0; t  < counter_triangles; t++  )
    for( int vi = 0; vi <                 3; vi++ )
    {
      int v = data_triangle_vertices[t][vi];
      
      assert( 0 <= v && v < counter_vertices );
      
      if( data_vertex_firstparent_triangle[v] == nullindex ) {
        
        data_vertex_firstparent_triangle[v] = t;
        
      } else {
        
        int old_first_parent = data_vertex_firstparent_triangle[v];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_triangles );
        assert( data_triangle_nextparents_of_vertices[ t ][ vi ] == nullindex );
        
        data_vertex_firstparent_triangle[v] = t;
        data_triangle_nextparents_of_vertices[ t ][ vi ] = old_first_parent;
        
      }
      
      assert( data_vertex_firstparent_triangle[v] != nullindex );
      assert( 0 <= data_vertex_firstparent_triangle[v] && data_vertex_firstparent_triangle[v] < counter_triangles );  
    }
    
    /* 4. For each vertex, set the first parent edge and the neighboring parent edges */
    
    for( int e =  0; e  < counter_edges; e++  )
    for( int vi = 0; vi <             2; vi++ )
    {
      int v = data_edge_vertices[e][vi];
      
      assert( 0 <= v && v < counter_vertices );
      
      if( data_vertex_firstparent_edge[v] == nullindex ) {
        
        data_vertex_firstparent_edge[v] = e;
        
      } else {
        
        int old_first_parent = data_vertex_firstparent_edge[v];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_edges );
        assert( data_edge_nextparents_of_vertices[ e ][ vi ] == nullindex );
        
        data_vertex_firstparent_edge[v] = e;
        data_edge_nextparents_of_vertices[ e ][ vi ] = old_first_parent;
        
      }
      
      assert( data_vertex_firstparent_edge[v] != nullindex );
      assert( 0 <= data_vertex_firstparent_edge[v] && data_vertex_firstparent_edge[v] < counter_edges );  
    }
    
    /* 5. For each edge, set the first parent triangle and the neighboring parent triangles */
    
    for( int t = 0; t < counter_triangles; t++ )
    for( int e = 0; e <     counter_edges; e++ )
    {
      int voe0 = data_edge_vertices[e][0];
      int voe1 = data_edge_vertices[e][1];
      
      int vot0 = data_triangle_vertices[t][0];
      int vot1 = data_triangle_vertices[t][1];
      int vot2 = data_triangle_vertices[t][2];
      
      int eot = nullindex;
      
      if( voe0 == vot0 && voe1 == vot1 ) eot = 0;
      if( voe0 == vot0 && voe1 == vot2 ) eot = 1;
      if( voe0 == vot1 && voe1 == vot2 ) eot = 2;
      
      if( eot == nullindex ) continue;
      
      data_triangle_edges[t][eot] = e;
      
      if( data_edge_firstparent_triangle[e] == nullindex ) {
        
        data_edge_firstparent_triangle[e] = t;
        
      } else {
        
        int old_first_parent = data_edge_firstparent_triangle[e];
        
        assert( 0 <= old_first_parent && old_first_parent < counter_triangles );
        
        data_edge_firstparent_triangle[e] = t;
        
        assert( data_triangle_nextparents_of_edges[ t ][ eot ] == nullindex );
        data_triangle_nextparents_of_edges[ t ][ eot ] = old_first_parent;
        
      }
      
      assert( data_edge_firstparent_triangle[e] != nullindex );
      assert( 0 <= data_edge_firstparent_triangle[e] && data_edge_firstparent_triangle[e] < counter_triangles );  
    }
    
    
    /* Coda: set the flags to the Null flag */
    flags_edges.resize   ( counter_edges,    SimplexFlag::SimplexFlagInvalid );
    flags_vertices.resize( counter_vertices, SimplexFlag::SimplexFlagInvalid );
    
    for( int t =  0; t  <  counter_triangles; t++  ) flags_triangles.at( t ) = SimplexFlag::SimplexFlagNull;
    for( int e =  0; e  <  counter_edges;     e++  ) flags_edges.at( e )     = SimplexFlag::SimplexFlagNull;
    for( int v =  0; v  <  counter_vertices;  v++  ) flags_vertices.at( v )  = SimplexFlag::SimplexFlagNull;
    
    
    MeshSimplicial2D::check();
}


MeshSimplicial2D::MeshSimplicial2D( 
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
)
:
    Mesh( 2, outerdim ),
    
    counter_triangles( triangle_vertices.size() ),
    counter_edges( edge_vertices.size() ),
    counter_vertices( vertex_firstparent_triangle.size() ),
    
    data_triangle_edges( triangle_edges ),
    data_edge_firstparent_triangle( edge_firstparent_triangle ),
    data_triangle_nextparents_of_edges( triangle_nextparents_of_edges ),
    
    data_triangle_vertices( triangle_vertices ),
    data_vertex_firstparent_triangle( vertex_firstparent_triangle ),
    data_triangle_nextparents_of_vertices( triangle_nextparents_of_vertices ),
    
    data_edge_vertices( edge_vertices ),
    data_vertex_firstparent_edge( vertex_firstparent_edge ),
    data_edge_nextparents_of_vertices( edge_nextparents_of_vertices ),
    
    flags_triangles( counter_triangles, SimplexFlag::SimplexFlagNull ),
    flags_edges    ( counter_edges,     SimplexFlag::SimplexFlagNull ),
    flags_vertices ( counter_vertices,  SimplexFlag::SimplexFlagNull )

{
    
    getcoordinates() = coords;
    
    /* set the flags to the Null flag */
    flags_edges.resize   ( counter_edges    );
    flags_vertices.resize( counter_vertices );
    
    for( int t =  0; t  <  counter_triangles; t++  ) flags_triangles.at( t ) = SimplexFlag::SimplexFlagNull;
    for( int e =  0; e  <  counter_edges;     e++  ) flags_edges.at( e )     = SimplexFlag::SimplexFlagNull;
    for( int v =  0; v  <  counter_vertices;  v++  ) flags_vertices.at( v )  = SimplexFlag::SimplexFlagNull;

    MeshSimplicial2D::check();
}


MeshSimplicial2D::~MeshSimplicial2D()
{
    MeshSimplicial2D::check();
}


bool MeshSimplicial2D::is_equal_to( const MeshSimplicial2D& mesh ) const 
{
  return counter_triangles == mesh.counter_triangles
         &&
         counter_edges == mesh.counter_edges
         &&
         counter_vertices == mesh.counter_vertices
         &&
         data_triangle_edges == mesh.data_triangle_edges
         &&
         data_edge_firstparent_triangle == mesh.data_edge_firstparent_triangle
         &&
         data_triangle_nextparents_of_edges == mesh.data_triangle_nextparents_of_edges
         &&
         data_triangle_vertices == mesh.data_triangle_vertices 
         &&
         data_vertex_firstparent_triangle == mesh.data_vertex_firstparent_triangle
         &&
         data_triangle_nextparents_of_vertices == mesh.data_triangle_nextparents_of_vertices
         &&
         data_edge_vertices == mesh.data_edge_vertices 
         &&
         data_vertex_firstparent_edge == mesh.data_vertex_firstparent_edge
         &&
         data_edge_nextparents_of_vertices == mesh.data_edge_nextparents_of_vertices
         &&
         getinnerdimension() == mesh.getinnerdimension()
         &&
         getouterdimension() == mesh.getouterdimension()
         &&
         getcoordinates() == mesh.getcoordinates()
         &&
         flags_triangles == mesh.flags_triangles
         &&
         flags_edges == mesh.flags_edges
         &&
         flags_vertices == mesh.flags_vertices
         &&
         true;
}





void MeshSimplicial2D::check() const
{
    
    #if defined(DO_NOT_CHECK_MESHES)
    #warning Disabled check for 2D Simplicial Mesh 
    return;
    #endif
    
    #ifdef NDEBUG
    return;
    #else
    
    /****************************/
    /* 1. Check the array sizes */
    /****************************/
    
    assert( counter_triangles == data_triangle_edges.size() );
    assert( counter_triangles == data_triangle_nextparents_of_edges.size() );
    assert( counter_triangles == data_triangle_vertices.size() );
    assert( counter_triangles == data_triangle_nextparents_of_vertices.size() );
    assert( counter_edges == data_edge_vertices.size() );
    assert( counter_edges == data_edge_nextparents_of_vertices.size() );
    assert( counter_edges == data_edge_firstparent_triangle.size() );
    assert( counter_vertices == data_vertex_firstparent_edge.size() );
    assert( counter_vertices == data_vertex_firstparent_triangle.size() );
    
    assert( count_vertices() == getcoordinates().getnumber() );
    
    assert( counter_triangles == flags_triangles.size() );
    assert( counter_edges     == flags_edges.size()     );
    assert( counter_vertices  == flags_vertices.size()  );

    
    
    /********************************************************************************/
    /* 2. check that the internal data of each simplex make sense on each dimension */
    /********************************************************************************/
    
    /* each triangle: each edge is a valid index */
    /* each triangle: each edge is unique        */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
        
        assert( data_triangle_edges[t][0] != nullindex );
        assert( data_triangle_edges[t][1] != nullindex );
        assert( data_triangle_edges[t][2] != nullindex );
        
        assert( 0 <= data_triangle_edges[t][0] && data_triangle_edges[t][0] < counter_edges );
        assert( 0 <= data_triangle_edges[t][1] && data_triangle_edges[t][1] < counter_edges );
        assert( 0 <= data_triangle_edges[t][2] && data_triangle_edges[t][2] < counter_edges );
        
        assert( data_triangle_edges[t][0] != data_triangle_edges[t][1] );
        assert( data_triangle_edges[t][0] != data_triangle_edges[t][2] );
        assert( data_triangle_edges[t][1] != data_triangle_edges[t][2] );
        
    }
    
    /* each triangle: each vertex is a valid index */
    /* each triangle: each vertex is unique        */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
        
        assert( data_triangle_vertices[t][0] != nullindex );
        assert( data_triangle_vertices[t][1] != nullindex );
        assert( data_triangle_vertices[t][2] != nullindex );
        
        assert( 0 <= data_triangle_vertices[t][0] && data_triangle_vertices[t][0] < counter_vertices );
        assert( 0 <= data_triangle_vertices[t][1] && data_triangle_vertices[t][1] < counter_vertices );
        assert( 0 <= data_triangle_vertices[t][2] && data_triangle_vertices[t][2] < counter_vertices );
        
        assert( data_triangle_vertices[t][0] != data_triangle_vertices[t][1] );
        assert( data_triangle_vertices[t][0] != data_triangle_vertices[t][2] );
        assert( data_triangle_vertices[t][1] != data_triangle_vertices[t][2] );
        
    }
    
    /* each edge: each vertex is a valid index */
    /* each edge: each vertex is unique        */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        
        assert( data_edge_vertices[e][0] != nullindex );
        assert( data_edge_vertices[e][1] != nullindex );
        assert( 0 <= data_edge_vertices[e][0] && data_edge_vertices[e][0] < counter_vertices );
        assert( 0 <= data_edge_vertices[e][1] && data_edge_vertices[e][0] < counter_vertices );
        assert( data_edge_vertices[e][0] != data_edge_vertices[e][1] );
        
    }
    
    
    
    
    /****************************************/
    /* 2. check the uniqueness of simplices */
    /****************************************/
    
    /* check that all edges are unique, even up to permutation */
    
    for( int e1 = 0; e1 < counter_edges; e1++ )
    for( int e2 = 0; e2 < counter_edges; e2++ )
    {
        if( e1 == e2 ) continue;
        
        assert( data_edge_vertices[e1][0] != data_edge_vertices[e2][0] || data_edge_vertices[e1][1] != data_edge_vertices[e2][1] );
        assert( data_edge_vertices[e1][0] != data_edge_vertices[e2][1] || data_edge_vertices[e1][1] != data_edge_vertices[e2][0] );
    }
    
    
    /* check that all triangles are unique, even up to permutation    */
    /* check both the vertices and the edges listed for each triangle */
    
    for( int t1 = 0; t1 < counter_triangles; t1++ )
    for( int t2 = 0; t2 < counter_triangles; t2++ )
    {
        if( t1 == t2 ) continue;
        
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][0] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][1] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][2] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][0] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][2] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][1] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][1] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][0] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][2] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][1] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][2] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][0] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][2] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][0] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][1] );
        assert( data_triangle_vertices[t1][0] != data_triangle_vertices[t2][2] || data_triangle_vertices[t1][1] != data_triangle_vertices[t2][1] || data_triangle_vertices[t1][2] != data_triangle_vertices[t2][0] );
        
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][0] || data_triangle_edges[t1][1] != data_triangle_edges[t2][1] || data_triangle_edges[t1][2] != data_triangle_edges[t2][2] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][0] || data_triangle_edges[t1][1] != data_triangle_edges[t2][2] || data_triangle_edges[t1][2] != data_triangle_edges[t2][1] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][1] || data_triangle_edges[t1][1] != data_triangle_edges[t2][0] || data_triangle_edges[t1][2] != data_triangle_edges[t2][2] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][1] || data_triangle_edges[t1][1] != data_triangle_edges[t2][2] || data_triangle_edges[t1][2] != data_triangle_edges[t2][0] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][2] || data_triangle_edges[t1][1] != data_triangle_edges[t2][0] || data_triangle_edges[t1][2] != data_triangle_edges[t2][1] );
        assert( data_triangle_edges[t1][0] != data_triangle_edges[t2][2] || data_triangle_edges[t1][1] != data_triangle_edges[t2][1] || data_triangle_edges[t1][2] != data_triangle_edges[t2][0] );
        
    }
    
    
    
    
    
    
    /********************************************************/
    /* 3. check the data of each simplex accross dimensions */
    /********************************************************/
    
    /* each triangle: each edge is listed correctly */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
        
        assert( data_edge_vertices[ data_triangle_edges[t][0] ][0] == data_triangle_vertices[t][0] );
        assert( data_edge_vertices[ data_triangle_edges[t][0] ][1] == data_triangle_vertices[t][1] );
        assert( data_edge_vertices[ data_triangle_edges[t][1] ][0] == data_triangle_vertices[t][0] );
        assert( data_edge_vertices[ data_triangle_edges[t][1] ][1] == data_triangle_vertices[t][2] );
        assert( data_edge_vertices[ data_triangle_edges[t][2] ][0] == data_triangle_vertices[t][1] );
        assert( data_edge_vertices[ data_triangle_edges[t][2] ][1] == data_triangle_vertices[t][2] );
        
    }
    
    
    
    
    
    
    /**************************************************************/
    /* 4. Check that the first parents are set and actual parents */
    /**************************************************************/
    
    /* each first parent triangle of an edge: first parent is non-null and a valid parent */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        int p = data_edge_firstparent_triangle[e];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_triangles );
        
        assert( data_triangle_edges[p][0] == e || data_triangle_edges[p][1] == e || data_triangle_edges[p][2] == e );
    }
    
    /* each first parent triangle of a vertex: first parent is non-null and a valid parent */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent_triangle[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_triangles );
        
        assert( data_triangle_vertices[p][0] == v || data_triangle_vertices[p][1] == v || data_triangle_vertices[p][2] == v );
    }
    
    /* each first parent edge of a vertex: first parent is non-null and a valid parent */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent_edge[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_edges );
        
        assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
    }
    
    
    
    /*****************************************************/
    /* 5. Check that the next parents are actual parents */
    /*****************************************************/
    
    /* each triangle: the next parents of edges are unique           */
    /* each triangle: the next parents of edges are actually parents */
    
    for( int t = 0; t < counter_triangles; t++ )
    for( int ei = 0; ei < 3; ei++ )
    {
        
        if( data_triangle_nextparents_of_edges[t][ei] != nullindex )
            assert( 0 <= data_triangle_nextparents_of_edges[t][ei] && data_triangle_nextparents_of_edges[t][ei] < counter_triangles );
    
        if( data_triangle_nextparents_of_edges[t][ei] != nullindex )
            assert( data_triangle_edges[ data_triangle_nextparents_of_edges[t][ei] ][0] == data_triangle_edges[t][ei] 
                    ||
                    data_triangle_edges[ data_triangle_nextparents_of_edges[t][ei] ][1] == data_triangle_edges[t][ei] 
                    ||
                    data_triangle_edges[ data_triangle_nextparents_of_edges[t][ei] ][2] == data_triangle_edges[t][ei] );
        
    }
    
    /* each triangle: the next parents of triangles are unique           */
    /* each triangle: the next parents of triangles are actually parents */
    
    for( int t = 0; t < counter_triangles; t++ )
    for( int vi = 0; vi < 3; vi++ )
    {
        
        if( data_triangle_nextparents_of_vertices[t][vi] != nullindex )
            assert( 0 <= data_triangle_nextparents_of_vertices[t][vi] && data_triangle_nextparents_of_vertices[t][vi] < counter_triangles );
    
        if( data_triangle_nextparents_of_vertices[t][vi] != nullindex )
            assert( data_triangle_vertices[ data_triangle_nextparents_of_vertices[t][vi] ][0] == data_triangle_vertices[t][vi] 
                    ||
                    data_triangle_vertices[ data_triangle_nextparents_of_vertices[t][vi] ][1] == data_triangle_vertices[t][vi] 
                    ||
                    data_triangle_vertices[ data_triangle_nextparents_of_vertices[t][vi] ][2] == data_triangle_vertices[t][vi] );
        
    }
    
    /* each edge: the next parents of vertices make sense           */
    /* each edge: the next parents of vertices are actually parents */
    
    for( int e = 0; e < counter_edges; e++ )
    for( int vi = 0; vi < 2; vi++ )
    {
        
        if( data_edge_nextparents_of_vertices[e][vi] != nullindex )
        assert( 0 <= data_edge_nextparents_of_vertices[e][vi] && data_edge_nextparents_of_vertices[e][vi] < counter_edges );
        
        if( data_edge_nextparents_of_vertices[e][vi] != nullindex )
        assert( data_edge_vertices[ data_edge_nextparents_of_vertices[e][vi] ][0] == data_edge_vertices[e][vi] 
                ||
                data_edge_vertices[ data_edge_nextparents_of_vertices[e][vi] ][1] == data_edge_vertices[e][vi] 
                );
        
    }

    
    
    
    
    /****************************************/
    /* 6. Check that all parents are listed */
    /****************************************/
    
    /* check that each parent triangle of an edge is listed as a parent */
    
    for( int t  = 0; t  < counter_triangles; t++ )
    for( int ei = 0; ei <                 3; ei++ )
    {
      
      int e = data_triangle_edges[t][ei];
      
      int p = data_edge_firstparent_triangle[e];
      
      assert( p != nullindex );
      
      while( p != t && p != nullindex )
        if( data_triangle_edges[p][0] == e )
          p = data_triangle_nextparents_of_edges[p][0];
        else if( data_triangle_edges[p][1] == e )
          p = data_triangle_nextparents_of_edges[p][1];
        else if( data_triangle_edges[p][2] == e )
          p = data_triangle_nextparents_of_edges[p][2];
        else
          unreachable();
        
      assert( p == t );
      
    }
    
    /* check that each parent triangle of a vertex is listed as a parent */
    
    for( int t  = 0; t  < counter_triangles; t++ )
    for( int vi = 0; vi <                 3; vi++ )
    {
      
      int v = data_triangle_vertices[t][vi];
      
      int p = data_vertex_firstparent_triangle[v];
      
      assert( p != nullindex );
      
      while( p != t && p != nullindex )
        if( data_triangle_vertices[p][0] == v )
          p = data_triangle_nextparents_of_vertices[p][0];
        else if( data_triangle_vertices[p][1] == v )
          p = data_triangle_nextparents_of_vertices[p][1];
        else if( data_triangle_vertices[p][2] == v )
          p = data_triangle_nextparents_of_vertices[p][2];
        else
          unreachable();
        
      assert( p == t );
      
    }
    
    /* check that each parent edge of a vertex is listed as a parent */
    
    for( int e  = 0; e  < counter_edges; e++ )
    for( int vi = 0; vi <             2; vi++ )
    {
      
      int v = data_edge_vertices[e][vi];
      
      int p = data_vertex_firstparent_edge[v];
      
      assert( p != nullindex );
      
      while( p != e && p != nullindex )
        p = data_edge_nextparents_of_vertices[p][ ( data_edge_vertices[p][0] == v ) ? 0 : 1 ];
      
      assert( p == e );
      
    }
    
    
    
    /*
     * check that all the flags are valid
     */
    
    for( int t  = 0; t  <  counter_triangles; t++ )
        assert( flags_triangles[t] != SimplexFlag::SimplexFlagInvalid );

    for( int e  = 0; e  <  counter_edges; e++ )
        assert( flags_edges[e] != SimplexFlag::SimplexFlagInvalid );

    for( int v  = 0; v  <  counter_vertices; v++ )
        assert( flags_vertices[v] != SimplexFlag::SimplexFlagInvalid );
    
    
    

    
    /*************************************/
    /* CHECK FROM BASE CLASS PERSPECTIVE */
    /*************************************/
    
    Mesh::check();
    #endif
}






std::string MeshSimplicial2D::text() const
{
    std::ostringstream os;
    
    os << "Triangulation of 2D Manifold!" << nl;
    
    os << counter_triangles << space << counter_edges << space << counter_vertices << nl;
    
    
    
    os << "Triangle edges" << nl;
    
    for( const auto& triple : data_triangle_edges )
      os << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    os << "Edge first parent triangles" << nl;
    
    for( int fp : data_edge_firstparent_triangle )
      os << fp << nl;
    
    os << "Triangle next parents of edges" << nl;
    
    for( const auto& triple : data_triangle_nextparents_of_edges )
      os << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    
    
    
    os << "Triangle vertices" << nl;
    
    for( const auto& triple : data_triangle_vertices )
      os << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    os << "Vertex first parent triangles" << nl;
    
    for( int fp : data_vertex_firstparent_triangle )
      os << fp << nl;
    
    os << "Triangle next parents of vertices" << nl;
    
    for( const auto& triple : data_triangle_nextparents_of_vertices )
      os << triple[0] << space << triple[1] << space << triple[2] << nl;
    
    
    
    os << "Edge vertices" << nl;
    
    for( const auto& duple : data_edge_vertices )
      os << duple[0] << space << duple[1] << nl;
    
    os << "Vertex first parents" << nl;
    
    for( int fp : data_vertex_firstparent_edge )
      os << fp << nl;
    
    os << "Edge next parents of vertices" << nl;
    
    for( const auto& duple : data_edge_nextparents_of_vertices )
      os << duple[0] << space << duple[1] << nl;
    
    os << "Finished printing" << nl;
    
    return os.str();
}






bool MeshSimplicial2D::has_dimension_counted( int dim ) const
{
    assert( 0 <= dim && dim <= 2 );
    return true;
}

int MeshSimplicial2D::count_simplices( int dim ) const
{
  if( dim == 0 )
    return count_vertices();
  else if( dim == 1 )
    return count_edges();
  else if( dim == 2 )
    return count_triangles();
  else
    unreachable();
}

bool MeshSimplicial2D::has_subsimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= 2 );
    return true;
}

IndexMap MeshSimplicial2D::getsubsimplices( int sup, int sub, int cell ) const
{
  assert( 0 <= sub && sub <= sup && sup <= 2 );
  assert( 0 <= cell );
  if( sup == 0 ) assert( cell <= count_vertices()  );
  if( sup == 1 ) assert( cell <= count_edges()     );
  if( sup == 2 ) assert( cell <= count_triangles() );
  
  if( sup == 2 && sub == 2 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_triangles()-1), { cell } );
    
  } else if( sup == 2 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    auto temp = get_triangle_edges(cell);
    return IndexMap( IndexRange(0,2), IndexRange( 0, count_edges()-1 ), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 2 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    auto temp = get_triangle_vertices(cell);
    return IndexMap( IndexRange(0,2), IndexRange( 0, count_vertices()-1 ) , std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 1 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_edges()-1), { cell } );
    
  } else if( sup == 1 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_edge_vertices(cell);
    return IndexMap( IndexRange(0,1), IndexRange( 0, count_vertices()-1 ) , std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 0 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_vertices()-1), { cell } );
    
  } else {
    
    unreachable();
    
  }
   
}

bool MeshSimplicial2D::has_supersimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= 2 );
    return true;
}

const std::vector<int> MeshSimplicial2D::getsupersimplices( int sup, int sub, int cell ) const
{
  
  if( sup == 2 && sub == 2 ) {
    
    assert( 0 <= cell && cell < count_triangles() );
    return { cell };
    
  } else if( sup == 2 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_triangle_parents_of_edge( cell ); 
    return temp; //return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 2 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_triangle_parents_of_vertex( cell ); 
    return temp; //return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 1 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    return { cell };
    
  } else if( sup == 1 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_edge_parents_of_vertex( cell ); 
    return temp; //return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 0 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    return { cell };
    
  } else {
    
    unreachable();
    
  }
  
}




SimplexFlag MeshSimplicial2D::get_flag( int dim, int cell ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    if( dim == 0 ) {
        assert( 0 <= cell && cell < count_vertices() );
        return flags_vertices[cell];
    } else if( dim == 1 ) {
        assert( 0 <= cell && cell < count_edges() );
        return flags_edges[cell];
    } else if( dim == 2 ) {
        assert( 0 <= cell && cell < count_triangles() );
        return flags_triangles[cell];
    } else {
        unreachable();
    }
}
        
void MeshSimplicial2D::set_flag( int dim, int cell, SimplexFlag flag )
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    if( dim == 0 ) {
        assert( 0 <= cell && cell < count_vertices() );
        flags_vertices[cell] = flag;
    } else if( dim == 1 ) {
        assert( 0 <= cell && cell < count_edges() );
        flags_edges[cell] = flag;
    } else if( dim == 2 ) {
        assert( 0 <= cell && cell < count_triangles() );
        flags_triangles[cell] = flag;
    } else {
        unreachable();
    }
}









/* Count number of elements */

int MeshSimplicial2D::count_triangles() const
{
    return counter_triangles;
}

int MeshSimplicial2D::count_edges() const
{
    return counter_edges;
}

int MeshSimplicial2D::count_vertices() const
{
    return counter_vertices;
}




/* subsimplex relation of triangles and edges */

bool MeshSimplicial2D::contains_triangle_edge( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    
    return ( data_triangle_edges[t][0] == e ) || ( data_triangle_edges[t][1] == e ) || ( data_triangle_edges[t][2] == e );
} 

int MeshSimplicial2D::indexof_triangle_edge( int t, int e ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= e && e < counter_edges );
    if     ( data_triangle_edges[t][0] == e ) return 0;
    else if( data_triangle_edges[t][1] == e ) return 1;
    else if( data_triangle_edges[t][2] == e ) return 2;
    else                                      unreachable();
} 

int MeshSimplicial2D::get_triangle_edge( int t, int ei ) const
{
    assert( 0 <= t  && t  < counter_triangles );
    assert( 0 <= ei && ei < 3 );
    return data_triangle_edges[t][ei];
} 

const std::array<int,3> MeshSimplicial2D::get_triangle_edges( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_edges[t];
} 



/* subsimplex relation of triangle and vertices */

bool MeshSimplicial2D::contains_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_triangle_vertices[t][0] == v ) || ( data_triangle_vertices[t][1] == v ) || ( data_triangle_vertices[t][2] == v );
} 

int MeshSimplicial2D::indexof_triangle_vertex( int t, int v ) const
{
    assert( 0 <= t && t < counter_triangles );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_triangle_vertices[t][0] == v ) return 0;
    else if( data_triangle_vertices[t][1] == v ) return 1;
    else if( data_triangle_vertices[t][2] == v ) return 2;
    else                                         unreachable();
} 

int MeshSimplicial2D::get_triangle_vertex( int t, int vi ) const
{
    assert( 0 <= t  && t  < counter_triangles );
    assert( 0 <= vi && vi < 3 );
    return data_triangle_vertices[t][vi];
} 

const std::array<int,3> MeshSimplicial2D::get_triangle_vertices( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    return data_triangle_vertices[t];
} 




/* subsimplex relation of edges and vertices */

bool MeshSimplicial2D::contains_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_edge_vertices[e][0] == v ) || ( data_edge_vertices[e][1] == v );
} 

int MeshSimplicial2D::indexof_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_edge_vertices[e][0] == v ) return 0;
    else if( data_edge_vertices[e][1] == v ) return 1;
    else                                     unreachable();
} 

int MeshSimplicial2D::get_edge_vertex( int e, int vi ) const
{
    assert( 0 <= e  && e  < counter_edges );
    assert( 0 <= vi && vi < 2 );
    return data_edge_vertices[e][vi];
} 

const std::array<int,2> MeshSimplicial2D::get_edge_vertices( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    return data_edge_vertices[e];
} 





/* triangle parents of a edge */

int MeshSimplicial2D::count_edge_triangle_parents( int e ) const
{
  return get_triangle_parents_of_edge( e ).size();
}

int MeshSimplicial2D::get_edge_firstparent_triangle( int e ) const
{
  assert( 0 <= e && e < counter_edges );
  return data_edge_firstparent_triangle[ e ];
}

int MeshSimplicial2D::get_edge_nextparent_triangle( int e, int t ) const
{
  assert( 0 <= e && e < counter_edges );
  assert( 0 <= t && t < counter_triangles );
  
  if( data_triangle_edges[t][0] == e )
    return data_triangle_nextparents_of_edges[t][0];
  else if( data_triangle_edges[t][1] == e )
    return data_triangle_nextparents_of_edges[t][1];
  else if( data_triangle_edges[t][2] == e )
    return data_triangle_nextparents_of_edges[t][2];
  else
    unreachable();
}

int MeshSimplicial2D::get_triangle_nextparent_of_edge( int t, int ei ) const
{
  assert( 0 <= t  && t  < counter_triangles );
  assert( 0 <= ei && ei < 3 );
  return data_triangle_nextparents_of_edges[t][ei];
}

bool MeshSimplicial2D::is_triangle_edge_parent( int t, int e ) const
{
  assert( 0 <= e && e < counter_edges );
  assert( 0 <= t && t < counter_triangles );
  return data_triangle_edges[t][0] == e || data_triangle_edges[t][1] == e || data_triangle_edges[t][2] == e;
}

int MeshSimplicial2D::indexof_triangle_edge_parent( int t, int e ) const
{
  assert( 0 <= e && e < count_edges() );
  std::vector<int> triangles = get_triangle_parents_of_edge( e );
  
  auto iter = std::find( triangles.begin(), triangles.end(), t ); 
  assert( iter != triangles.end() );
  
  return iter - triangles.begin();
}

std::vector<int> MeshSimplicial2D::get_triangle_parents_of_edge( int e ) const
{
  assert( 0 <= e && e < count_edges() );
  
  std::vector<int> ret;
  
  int t = get_edge_firstparent_triangle(e);
  while( t != nullindex )
  {
    ret.push_back(t);
    t = this->get_edge_nextparent_triangle(e,t);
  }
  
  return ret;
}


/* triangle parents of a vertex */

int MeshSimplicial2D::count_vertex_triangle_parents( int v ) const
{
  return get_triangle_parents_of_vertex( v ).size();
}

int MeshSimplicial2D::get_vertex_firstparent_triangle( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent_triangle[ v ];
}

int MeshSimplicial2D::get_vertex_nextparent_triangle( int v, int t ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= t && t < counter_triangles );
  
  if( data_triangle_vertices[t][0] == v )
    return data_triangle_nextparents_of_vertices[t][0];
  else if( data_triangle_vertices[t][1] == v )
    return data_triangle_nextparents_of_vertices[t][1];
  else if( data_triangle_vertices[t][2] == v )
    return data_triangle_nextparents_of_vertices[t][2];
  else
    unreachable();
}

int MeshSimplicial2D::get_triangle_nextparent_of_vertex( int t, int vi ) const
{
  assert( 0 <= t  && t  < counter_triangles );
  assert( 0 <= vi && vi < 3 );
  return data_triangle_nextparents_of_vertices[t][vi];
}

bool MeshSimplicial2D::is_triangle_vertex_parent( int t, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= t && t < counter_triangles );
  return data_triangle_vertices[t][0] == v || data_triangle_vertices[t][1] == v || data_triangle_vertices[t][2] == v;
}

int MeshSimplicial2D::indexof_triangle_vertex_parent( int t, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> triangles = get_triangle_parents_of_vertex( v );
  
  auto iter = std::find( triangles.begin(), triangles.end(), t ); 
  assert( iter != triangles.end() );
  
  return iter - triangles.begin();
}

std::vector<int> MeshSimplicial2D::get_triangle_parents_of_vertex( int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  
  std::vector<int> ret;
  
  int t = get_vertex_firstparent_triangle(v);
  while( t != nullindex )
  {
    ret.push_back(t);
    t = this->get_vertex_nextparent_triangle(v,t);
  }

  return ret;
}




/* edge parents of a vertex */

int MeshSimplicial2D::count_vertex_edge_parents( int v ) const
{
  return get_edge_parents_of_vertex( v ).size();
}

int MeshSimplicial2D::get_vertex_firstparent_edge( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent_edge[ v ];
}

int MeshSimplicial2D::get_vertex_nextparent_edge( int v, int e ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= e && e < counter_edges    );
  
  if( data_edge_vertices[e][0] == v )
    return data_edge_nextparents_of_vertices[e][0];
  else if( data_edge_vertices[e][1] == v )
    return data_edge_nextparents_of_vertices[e][1];
  else
    unreachable();
}

int MeshSimplicial2D::get_edge_nextparent_of_vertex( int e, int vi ) const
{
  assert( 0 <= e  && e  < counter_edges );
  assert( 0 <= vi && vi < 2 );
  return data_edge_nextparents_of_vertices[e][vi];
}

bool MeshSimplicial2D::is_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= e && e < counter_edges    );
  return data_edge_vertices[e][0] == v || data_edge_vertices[e][1] == v;
}

int MeshSimplicial2D::indexof_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> edges = get_edge_parents_of_vertex( v );
  
  auto iter = std::find( edges.begin(), edges.end(), e ); 
  assert( iter != edges.end() );
  
  return iter - edges.begin();
}

std::vector<int> MeshSimplicial2D::get_edge_parents_of_vertex( int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  
  std::vector<int> ret;
  
  int e = get_vertex_firstparent_edge(v);
  while( e != nullindex )
  {
    ret.push_back(e);
    e = this->get_vertex_nextparent_edge(v,e);
  }

  return ret;
}




/*
 * * * * * NEWEST VERTEX BISECTION
 */
    
void MeshSimplicial2D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    // check();
    
    /* 
     * DESCRIPTION
     * 
     * The Bisection works as follows:
     * 
     * triangle -> edge 
     * 
     *          T2E: 
     *          
     *            just fill in the data. 
     *          
     *          next parents + first parent:
     *            
     *            for outer edges:
     *            start with the first parent pointer and run over all parents,
     *            replacing the old triangle with the new triangle if need be.
     *            
     *            for inner edges: 
     *            fill in the data
     *            
     *            Finally, create parent data for the two new edges 
     *          
     * 
     * triangle -> vertex 
     *          
     *          T2V: 
     *          
     *            just fill in the data.
     *          
     *          next parents + first parent:
     *          
     *            for front,back,opposing vertex:
     *            start with the first parent pointer and run over all parents,
     *            replacing the old triangle with the new triangle if need be.
     *            
     *            Finally, create new data for new vertex.
     *          
     *  
     * edge -> vertex 
     * 
     *          E2V: 
     *          
     *            just fill in the data.
     *            
     *          next parents + first parent:
     *          
     *            for back vertex:
     *                  [nothing changes]
     *            for front vertex:
     *                  run over all parents and replace the old edge.
     *            for opposing vertex:
     *                  add the additional edge as first parent
     * 
     * 
     * */
    
    
    
    
    /*
     * DATA COLLECTION 
     */
    
    int e_back_vertex  = data_edge_vertices[ e ][ 0 ];
    int e_front_vertex = data_edge_vertices[ e ][ 1 ];
    
    std::vector<int> old_triangles = get_triangle_parents_of_edge( e );
    
    std::vector<int> localindex_of_refinementedge( old_triangles.size() );
    
    
    FloatVector midcoordinate = get_edge_midpoint( e );
    
    
    SimplexFlag e_flag = flags_edges[e];
    
    
    /*
     * ALLOCATE MEMORY FOR THE DATA  
     */
    
    data_triangle_nextparents_of_edges.resize( counter_triangles + old_triangles.size()    , { nullindex, nullindex, nullindex } );
    data_triangle_edges.resize               ( counter_triangles + old_triangles.size()    , { nullindex, nullindex, nullindex } );
    data_edge_firstparent_triangle.resize    ( counter_edges     + old_triangles.size() + 1,                           nullindex );
    
    data_triangle_nextparents_of_vertices.resize( counter_triangles + old_triangles.size(), { nullindex, nullindex, nullindex } );
    data_triangle_vertices.resize               ( counter_triangles + old_triangles.size(), { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_triangle.resize     ( counter_vertices  + 1,                                              nullindex );
    
    data_edge_nextparents_of_vertices.resize( counter_edges    + old_triangles.size() + 1, { nullindex, nullindex } );
    data_edge_vertices.resize               ( counter_edges    + old_triangles.size() + 1, { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1,                                       nullindex );
    
    
    flags_triangles.resize( counter_triangles + old_triangles.size(),     SimplexFlag::SimplexFlagInvalid );
    flags_edges.resize    ( counter_edges     + old_triangles.size() + 1, SimplexFlag::SimplexFlagInvalid );
    flags_vertices.resize ( counter_vertices  + 1,                        SimplexFlag::SimplexFlagInvalid );
    

    
    
    /*
     * FILL IN DATA
     */
    
    /* vertices of the bisected edge */
    data_edge_vertices[ e ][ 0 ] = e_back_vertex;
    data_edge_vertices[ e ][ 1 ] = counter_vertices;
      
    data_edge_vertices[ counter_edges ][ 0 ] = counter_vertices;
    data_edge_vertices[ counter_edges ][ 1 ] = e_front_vertex;
    
    /* next parent of back vertex stays the same */
    /* next parent of front vertex */
    data_edge_nextparents_of_vertices[ counter_edges ][ 1 ] = data_edge_nextparents_of_vertices[ e ][ 1 ];
    
    // edge parent list of back vertex stays the same 
    // run over the front vertex edge parent list and replace 'e' by 'counter_edges'
    if( data_vertex_firstparent_edge[ e_front_vertex ] == e ) {
      
      data_vertex_firstparent_edge[ e_front_vertex ] = counter_edges;
      
    } else {
      
      int current_edge = data_vertex_firstparent_edge[ e_front_vertex ];
      
      while( data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] != e )
        current_edge = data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ];
      
      data_edge_nextparents_of_vertices[ current_edge ][ indexof_edge_vertex( current_edge, e_front_vertex ) ] = counter_edges;
      
    }
    
    /* first and next parents of new vertex */
    data_vertex_firstparent_edge[ counter_vertices ] = e;
    data_edge_nextparents_of_vertices[ e ][ 1 ] = counter_edges;
    data_edge_nextparents_of_vertices[ counter_edges ][ 0 ] = nullindex;
    
    // no parent triangles yet for the new vertex; will be filled in below */
    data_vertex_firstparent_triangle[ counter_vertices ] = nullindex;
    
    // no parent triangles yet for the front edge; will be filled in below */
    data_edge_firstparent_triangle[ counter_edges ] = nullindex;
    
    
    /*
     * flags of the bisected edge, its immediate children, and the new vertex 
     */
    SimplexFlag flag_oldedge = flags_edges[e];
    
    flags_edges[ e             ] = flag_oldedge;
    flags_edges[ counter_edges ] = flag_oldedge;
    
    flags_vertices[ counter_vertices ] = flag_oldedge;
    
    
    
    for( int ot = 0; ot < old_triangles.size(); ot++ ) {
      
      int t_old = old_triangles[ ot ];
      int t_new = counter_triangles + ot;
      
      int t_e0 = data_triangle_edges[ t_old ][ 0 ];
      int t_e1 = data_triangle_edges[ t_old ][ 1 ];
      int t_e2 = data_triangle_edges[ t_old ][ 2 ];
      
      int t_v0 = data_triangle_vertices[ t_old ][ 0 ];
      int t_v1 = data_triangle_vertices[ t_old ][ 1 ];
      int t_v2 = data_triangle_vertices[ t_old ][ 2 ];
      
      int t_e_n0 = data_triangle_nextparents_of_edges[ t_old ][ 0 ];
      int t_e_n1 = data_triangle_nextparents_of_edges[ t_old ][ 1 ];
      int t_e_n2 = data_triangle_nextparents_of_edges[ t_old ][ 2 ];
      
      int t_v_n0 = data_triangle_nextparents_of_vertices[ t_old ][ 0 ];
      int t_v_n1 = data_triangle_nextparents_of_vertices[ t_old ][ 1 ];
      int t_v_n2 = data_triangle_nextparents_of_vertices[ t_old ][ 2 ];
      
      localindex_of_refinementedge[ ot ] = ( t_e0 == e ) ? 0 : ( ( t_e1 == e ) ? 1 : 2 );
      assert( data_triangle_edges[ t_old ][ localindex_of_refinementedge[ ot ] ] == e );
      
      
      /*
       * flags of old triangle
       */
      
      SimplexFlag flag_ot = flags_triangles[ old_triangles[ ot ] ];
      flags_triangles[ t_old ] = flag_ot;
      flags_triangles[ t_new ] = flag_ot;
      
      flags_edges[ counter_edges + 1 + ot ] = flag_ot;
      
      
      
      
      if( localindex_of_refinementedge[ ot ] == 0 ) { // 0 1 
        
        assert( t_v0 == e_back_vertex && t_v1 == e_front_vertex && e == t_e0 );
        
        /* triangle vertices */
        data_triangle_vertices[ t_old ][0] = t_v0;
        data_triangle_vertices[ t_old ][1] = counter_vertices;
        data_triangle_vertices[ t_old ][2] = t_v2;
        
        data_triangle_vertices[ t_new ][0] = counter_vertices;
        data_triangle_vertices[ t_new ][1] = t_v1;
        data_triangle_vertices[ t_new ][2] = t_v2;
        
        /* triangle edges */
        data_triangle_edges[ t_old ][0] = e;
        data_triangle_edges[ t_old ][1] = t_e1;
        data_triangle_edges[ t_old ][2] = counter_edges + 1 + ot;
        
        data_triangle_edges[ t_new ][0] = counter_edges;
        data_triangle_edges[ t_new ][1] = counter_edges + 1 + ot;
        data_triangle_edges[ t_new ][2] = t_e2;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + ot ][0] = counter_vertices;
        data_edge_vertices[ counter_edges + 1 + ot ][1] = t_v2;
        
        
        /* triangle vertex neighbors */
        data_triangle_nextparents_of_vertices[ t_old ][0] = t_v_n0;
        data_triangle_nextparents_of_vertices[ t_old ][1] = counter_triangles + ot;
        data_triangle_nextparents_of_vertices[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_vertices[ t_new ][0] = data_vertex_firstparent_triangle[ counter_vertices ];
        data_triangle_nextparents_of_vertices[ t_new ][1] = t_v_n1;
        data_triangle_nextparents_of_vertices[ t_new ][2] = t_v_n2;
        
        data_vertex_firstparent_triangle[ counter_vertices ] = t_old;
        
        /* triangle edge neighbors */
        data_triangle_nextparents_of_edges[ t_old ][0] = t_e_n0;
        data_triangle_nextparents_of_edges[ t_old ][1] = t_e_n1;
        data_triangle_nextparents_of_edges[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_edges[ t_new ][0] = data_edge_firstparent_triangle[ counter_edges ];
        data_triangle_nextparents_of_edges[ t_new ][1] = nullindex;
        data_triangle_nextparents_of_edges[ t_new ][2] = t_e_n2;
        
        data_edge_firstparent_triangle[ counter_edges ] = t_new;
        
        data_edge_firstparent_triangle[ counter_edges + 1 + ot ] = t_old;
        
        /* run over the front outer edge parent triangle list and replace 't_old' by 't_new' */
        if( data_edge_firstparent_triangle[ t_e2 ] == t_old ) {
          data_edge_firstparent_triangle[ t_e2 ] = t_new;
        } else {
          int current_triangle = data_edge_firstparent_triangle[ t_e2 ];
          while( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != t_old 
                 &&
                 data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex )
            current_triangle = data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ];
          assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex );
          data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] = t_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][0] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ t_v2 ];
        data_vertex_firstparent_edge[ t_v2 ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][1] = opposing_vertex_firstparent_edge;
        
        
        
      } else if( localindex_of_refinementedge[ ot ] == 1 ) { // 0 2 
        
        assert( t_v0 == e_back_vertex && t_v2 == e_front_vertex && e == t_e1 );
        
        /* triangle vertices */
        data_triangle_vertices[ t_old ][0] = t_v0;
        data_triangle_vertices[ t_old ][1] = t_v1;
        data_triangle_vertices[ t_old ][2] = counter_vertices;
        
        data_triangle_vertices[ t_new ][0] = t_v1;
        data_triangle_vertices[ t_new ][1] = counter_vertices;
        data_triangle_vertices[ t_new ][2] = t_v2;
        
        /* triangle edges */
        data_triangle_edges[ t_old ][0] = t_e0;
        data_triangle_edges[ t_old ][1] = e;
        data_triangle_edges[ t_old ][2] = counter_edges + 1 + ot;
        
        data_triangle_edges[ t_new ][0] = counter_edges + 1 + ot;
        data_triangle_edges[ t_new ][1] = t_e2;
        data_triangle_edges[ t_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + ot ][0] = t_v1;
        data_edge_vertices[ counter_edges + 1 + ot ][1] = counter_vertices;
        
        
        /* triangle vertex neighbors */
        data_triangle_nextparents_of_vertices[ t_old ][0] = t_v_n0;
        data_triangle_nextparents_of_vertices[ t_old ][1] = counter_triangles + ot;
        data_triangle_nextparents_of_vertices[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_vertices[ t_new ][0] = t_v_n1;
        data_triangle_nextparents_of_vertices[ t_new ][1] = data_vertex_firstparent_triangle[ counter_vertices ];
        data_triangle_nextparents_of_vertices[ t_new ][2] = t_v_n2;
        
        data_vertex_firstparent_triangle[ counter_vertices ] = t_old;
        
        /* triangle edge neighbors */
        data_triangle_nextparents_of_edges[ t_old ][0] = t_e_n0;
        data_triangle_nextparents_of_edges[ t_old ][1] = t_e_n1;
        data_triangle_nextparents_of_edges[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_edges[ t_new ][0] = nullindex;
        data_triangle_nextparents_of_edges[ t_new ][1] = t_e_n2;
        data_triangle_nextparents_of_edges[ t_new ][2] = data_edge_firstparent_triangle[ counter_edges ];
        
        data_edge_firstparent_triangle[ counter_edges ] = t_new;
        
        data_edge_firstparent_triangle[ counter_edges + 1 + ot ] = t_old;
        
        /* run over the front outer edge parent triangle list and replace 't_old' by 't_new' */
        if( data_edge_firstparent_triangle[ t_e2 ] == t_old ) {
          data_edge_firstparent_triangle[ t_e2 ] = t_new;
        } else {
          int current_triangle = data_edge_firstparent_triangle[ t_e2 ];
          while( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != t_old 
                 &&
                 data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex )
            current_triangle = data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ];
          assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex );
          data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] = t_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][1] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ t_v1 ];
        data_vertex_firstparent_edge[ t_v1 ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][0] = opposing_vertex_firstparent_edge;
        
        
      } else if( localindex_of_refinementedge[ ot ] == 2 ) { // 1 2 
        
        assert( t_v1 == e_back_vertex && t_v2 == e_front_vertex && e == t_e2 );
        
        /* triangle vertices */
        data_triangle_vertices[ t_old ][0] = t_v0;
        data_triangle_vertices[ t_old ][1] = t_v1;
        data_triangle_vertices[ t_old ][2] = counter_vertices;
        
        data_triangle_vertices[ t_new ][0] = t_v0;
        data_triangle_vertices[ t_new ][1] = counter_vertices;
        data_triangle_vertices[ t_new ][2] = t_v2;
        
        /* triangle edges */
        data_triangle_edges[ t_old ][0] = t_e0;
        data_triangle_edges[ t_old ][1] = counter_edges + 1 + ot;
        data_triangle_edges[ t_old ][2] = e;
        
        data_triangle_edges[ t_new ][0] = counter_edges + 1 + ot;
        data_triangle_edges[ t_new ][1] = t_e1;
        data_triangle_edges[ t_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + ot ][0] = t_v0;
        data_edge_vertices[ counter_edges + 1 + ot ][1] = counter_vertices;
        
        
        /* triangle vertex neighbors */
        data_triangle_nextparents_of_vertices[ t_old ][0] = counter_triangles + ot;
        data_triangle_nextparents_of_vertices[ t_old ][1] = t_v_n1;
        data_triangle_nextparents_of_vertices[ t_old ][2] = counter_triangles + ot;
        
        data_triangle_nextparents_of_vertices[ t_new ][0] = t_v_n0;
        data_triangle_nextparents_of_vertices[ t_new ][1] = data_vertex_firstparent_triangle[ counter_vertices ];
        data_triangle_nextparents_of_vertices[ t_new ][2] = t_v_n2;
        
        data_vertex_firstparent_triangle[ counter_vertices ] = t_old;
        
        /* triangle edge neighbors */
        data_triangle_nextparents_of_edges[ t_old ][0] = t_e_n0;
        data_triangle_nextparents_of_edges[ t_old ][1] = counter_triangles + ot;
        data_triangle_nextparents_of_edges[ t_old ][2] = t_e_n2;
        
        data_triangle_nextparents_of_edges[ t_new ][0] = nullindex;
        data_triangle_nextparents_of_edges[ t_new ][1] = t_e_n1;
        data_triangle_nextparents_of_edges[ t_new ][2] = data_edge_firstparent_triangle[ counter_edges ];
        
        data_edge_firstparent_triangle[ counter_edges ] = t_new;
        
        data_edge_firstparent_triangle[ counter_edges + 1 + ot ] = t_old;
        
        /* run over the front outer edge parent triangle list and replace 't_old' by 't_new' */
        if( data_edge_firstparent_triangle[ t_e1 ] == t_old ) {
          data_edge_firstparent_triangle[ t_e1 ] = t_new;
        } else {
          int current_triangle = data_edge_firstparent_triangle[ t_e1 ];
          while( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != t_old 
                 &&
                 data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != nullindex )
            current_triangle = data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ];
          assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != nullindex );
          data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] = t_new;
        }
        
        
        /* add new parent edge for the new vertex */
        int new_vertex_firstparent_edge = data_vertex_firstparent_edge[ counter_vertices ];
        data_vertex_firstparent_edge[ counter_vertices ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][1] = new_vertex_firstparent_edge;
        
        /* add new parent edge for the opposing vertex */
        int opposing_vertex_firstparent_edge = data_vertex_firstparent_edge[ t_v0 ];
        data_vertex_firstparent_edge[ t_v0 ] = counter_edges + 1 + ot;
        data_edge_nextparents_of_vertices[ counter_edges + 1 + ot ][0] = opposing_vertex_firstparent_edge;
        
      } else {
        
        unreachable();
        
      } 
      
    }
    
    
    /* Run over the front vertex' parent triangles and conduct manipulations */

    int* pointer_to_index = &data_vertex_firstparent_triangle[ e_front_vertex ];
    
    while( *pointer_to_index != nullindex ) { 
      
      std::vector<int>::iterator it = std::find( old_triangles.begin(), old_triangles.end(), *pointer_to_index );
      
      if( it != old_triangles.end() ) {
        
        assert( *pointer_to_index == *it );
        assert( *it == old_triangles[ it - old_triangles.begin() ] );
        
        *pointer_to_index = counter_triangles + ( it - old_triangles.begin() );
        
      }
      
      int localindex_of_front_vertex = nullindex;
      if( data_triangle_vertices[ *pointer_to_index ][ 0 ] == e_front_vertex ) localindex_of_front_vertex = 0;
      if( data_triangle_vertices[ *pointer_to_index ][ 1 ] == e_front_vertex ) localindex_of_front_vertex = 1;
      if( data_triangle_vertices[ *pointer_to_index ][ 2 ] == e_front_vertex ) localindex_of_front_vertex = 2;
      assert( localindex_of_front_vertex != nullindex );
      
      pointer_to_index = &( data_triangle_nextparents_of_vertices[ *pointer_to_index ][ localindex_of_front_vertex ] );
      
    }
    
    getcoordinates().append( midcoordinate );
    
    
    /*
     *  UPDATE COUNTERS 
     */
    
    counter_triangles += old_triangles.size();
    counter_edges     += 1 + old_triangles.size();
    counter_vertices  += 1;
    
    
    
    
    
    /* Done */
    
    // check();
    
}






/*
 * * * * * LONGEST EDGE VERTEX BISECTION
 */

void MeshSimplicial2D::longest_edge_bisection_recursive( const std::vector<int>& edges )
{
    check();

    // 0. check the input
    
    for( const int& e : edges )
        assert( 0 <= e && e < counter_edges );


    // 1. run over the edges and apply the recursion

    const int old_vertex_count = counter_vertices;

    for( const int& e: edges ) {

        if( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count )
            continue;

        longest_edge_bisection_recursive( e );

    }

    // 3. check the result 
    
    for( const int& e : edges )
        assert( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count );
    
    check();
}


void MeshSimplicial2D::longest_edge_bisection_recursive( const int e )
{
    assert( 0 <= e && e < counter_edges );

    int longest_edge = nullindex;
        
    do{

        longest_edge = e;
        
        for( int t = get_edge_firstparent_triangle( e ); t != nullindex; t = get_edge_nextparent_triangle( e, t ) )
        for( int ei = 0; ei < 3; ei++ )
        {
            int other_e = get_triangle_edge( t, ei );
            
            if( e != other_e && get_edge_length( other_e ) > get_edge_length( e ) )
                longest_edge = other_e;
        }

        if( longest_edge != e )
            longest_edge_bisection_recursive( longest_edge );

    } while ( longest_edge != e );

    bisect_edge( e );

}



// void MeshSimplicial2D::longest_edge_bisection( std::vector<int> edges )
// {
//     check();
    
//     const int old_vertex_count = counter_vertices;
    
    
//     /* 0. check the input */
    
//     for( int& e : edges )
//         assert( 0 <= e && e < counter_edges );
    
    
//     /* 1. create stack for the edges to be bisected, and fill in first batch */
    
//     std::stack<int> todostack;
    
//     for( int& e : edges )
//         todostack.push( e );
    
// //     { /* remove duplicates */
// //         std::sort( todostack.begin(), todostack.end() );
// //         auto last = std::unique( todostack.begin(), todostack.end() );
// //         todostack.erase( last, todostack.end() );
// //     }
        
    
//     /* 2. conduct the main loop of the refinement algorithm */
    
//     while( ! todostack.empty() )
//     {
        
//         // as long the stack is not empty,
//         // pick the top edge and make the following case distinction
//         // a) the edge index belongs to an edge already bisected and can be ignored 
//         // b) if the top edge is bisection edge to all its neighbors, bisect and pop
//         // c) else, push the necessary edge of each parent simplex
        
//         int e = todostack.top();
        
//         // to check whether e belongs to an edge that has already been bisected,
//         // we check whether one of the vertices belongs to the new vertices 
        
//         if( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count ) {
            
//             todostack.pop();
        
//         } else {
            
//             // run over neighbor triangles and check for necessary edges 
            
//             for( int t = get_edge_firstparent_triangle( e ); t != nullindex; t = get_edge_nextparent_triangle( e, t ) )
//             {
                
//                 int e0 = get_triangle_edge( t, 0 ); Float l0 = get_edge_length( e0 );
//                 int e1 = get_triangle_edge( t, 1 ); Float l1 = get_edge_length( e1 );
//                 int e2 = get_triangle_edge( t, 2 ); Float l2 = get_edge_length( e2 );
                
//                 assert( e == e0 || e == e1 || e == e2 );
                
//                 int ne = nullindex;
//                 if( l0 >= l1 && l0 >= l2 ) ne = e0;
//                 if( l1 >= l0 && l1 >= l2 ) ne = e1;
//                 if( l2 >= l0 && l2 >= l1 ) ne = e2;
//                 assert( ne != nullindex );
                
//                 if( ne != e )
//                     todostack.push( ne );
                
//             }
            
//             // if top edge is still the same, bisect
            
//             if( e == todostack.top() )
//             {
//                 todostack.pop();
//                 bisect_edge( e );
//             }
            
//         }
        
//     }
    
    
//     /*
//     while( ! todostack.empty() )
//     {
        
//         LOG << todostack.back() << space << todostack.size();
        
//         // as long the stack is not empty,
//         // pick the top edge and make the following case distinction
//         // a) if the top edge is longer than its neighbors, bisect and pop
//         // b) else, push the longest edge of each parent simplex
        
//         int e = todostack.back();
    
//         Float length_e = get_edge_length( e );
            
//         // run over neighbor triangles and check for longer edges 
        
//         bool compatibly_divisible = true;
        
//         for( int t = get_edge_firstparent_triangle( e ); t != nullindex; t = get_edge_nextparent_triangle( e, t ) )
//         for( int ei = 0; ei < 3; ei++ )
//         {
//             int other_e = get_triangle_edge( t, ei );
            
//             if( e != other_e && get_edge_length( other_e ) > length_e ) {
//                 compatibly_divisible = false;
//                 todostack.push_back( other_e );
//             }
//         }
            
        
//         // if top edge is still the same, bisect 
    
//         if( compatibly_divisible )
//         {
//             assert( e == todostack.back() );
//             todostack.remove( e );
//             bisect_edge( e );
//             assert( std::find( todostack.begin(), todostack.end(), e ) == todostack.end() );
//         } else 
//             assert( e != todostack.back() );
        
    
//     }
//     */
    
//     // fiinished!
    
//     check();
// }







void MeshSimplicial2D::newest_vertex_bisection_recursive( const std::vector<int>& edges )
{
    check();

    // 0. check the input
    
    for( const int& e : edges )
        assert( 0 <= e && e < counter_edges );


    // 1. run over the edges and apply the recursion

    const int old_vertex_count = counter_vertices;

    for( const int& e: edges ) {

        if( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count )
            continue;

        newest_vertex_bisection_recursive( e );

    }

    // 3. check the result 
    
    for( const int& e : edges )
        assert( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count );
    
    check();
}


void MeshSimplicial2D::newest_vertex_bisection_recursive( const int e )
{
    assert( 0 <= e && e < counter_edges );

    int bisection_edge = nullindex;
        
    do{

        bisection_edge = nullindex;
        
        for( int t = get_edge_firstparent_triangle( e ); t != nullindex; t = get_edge_nextparent_triangle( e, t ) )
        {
            int v0 = get_triangle_vertex( t, 0 );
            int v1 = get_triangle_vertex( t, 1 );
            int v2 = get_triangle_vertex( t, 2 );
            int e0 = get_triangle_edge( t, 0 ); // 0 1
            int e1 = get_triangle_edge( t, 1 ); // 0 2
            int e2 = get_triangle_edge( t, 2 ); // 1 2 
            
            assert( e == e0 || e == e1 || e == e2 );
            
            if( v0 > v1 && v0 > v2 ) bisection_edge = e2;
            if( v1 > v0 && v1 > v2 ) bisection_edge = e1;
            if( v2 > v0 && v2 > v1 ) bisection_edge = e0;
            assert( bisection_edge != nullindex );
            
            if( bisection_edge != e ) break;       
        }
        
        if( bisection_edge != e )
            newest_vertex_bisection_recursive( bisection_edge );

    } while ( bisection_edge != e );

    bisect_edge( e );

}



// void MeshSimplicial2D::newest_vertex_bisection( std::vector<int> edges )
// {
//     check();
    
//     const int old_vertex_count = counter_vertices;
    
    
//     /* 0. check the input */
    
//     for( int& e : edges )
//         assert( 0 <= e && e < counter_edges );
    
    
//     /* 1. create stack for the edges to be bisected, and fill in first batch */
    
//     std::stack<int> todostack;
    
//     for( int& e : edges )
//         todostack.push( e ); // put e on top either by inserting or pulling it up!
        
    
//     /* 2. conduct the main loop of the refinement algorithm */
    
//     while( ! todostack.empty() )
//     {
        
//         // as long the stack is not empty,
//         // pick the top edge and make the following case distinction
//         // a) the edge index belongs to an edge already bisected and can be ignored 
//         // b) if the top edge is bisection edge to all its neighbors, bisect and pop
//         // c) else, push the necessary edge of each parent simplex
        
//         int e = todostack.top();
        
//         // to check whether e belongs to an edge that has already been bisected,
//         // we check whether one of the vertices belongs to the new vertices 
        
//         if( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count ) {
            
//             todostack.pop();
        
//         } else {
            
//             // run over neighbor triangles and check for necessary edges 
            
//             for(
//                 int t = get_edge_firstparent_triangle( e );
//                 t != nullindex; 
//                 t = get_edge_nextparent_triangle( e, t )
//             ) {
                
//                 int v0 = get_triangle_vertex( t, 0 );
//                 int v1 = get_triangle_vertex( t, 1 );
//                 int v2 = get_triangle_vertex( t, 2 );
//                 int e0 = get_triangle_edge( t, 0 ); // 0 1
//                 int e1 = get_triangle_edge( t, 1 ); // 0 2
//                 int e2 = get_triangle_edge( t, 2 ); // 1 2 
                
//                 assert( e == e0 || e == e1 || e == e2 );
                
//                 int ne = nullindex;
//                 if( v0 > v1 && v0 > v2 ) ne = e2;
//                 if( v1 > v0 && v1 > v2 ) ne = e1;
//                 if( v2 > v0 && v2 > v1 ) ne = e0;
//                 assert( ne != nullindex );
                
//                 if( ne != e )
//                     todostack.push( ne );
                
//             }
            
//             // if top edge is still the same, bisect
            
//             if( e == todostack.top() )
//             {
//                 todostack.pop();
//                 bisect_edge( e );
//             }
            
//         }
        
//     }
    
    
//     // fiinished!
    
//     check();
// }








void MeshSimplicial2D::uniformrefinement( int levels )
{
  assert( levels >= 0 );
  for( int l = 0; l < levels; l++ )
    uniformrefinement();
}



void MeshSimplicial2D::uniformrefinement()
{
    check();
    
    /* resize the arrays */
    
    data_triangle_edges.resize               ( 4 * counter_triangles                    , { nullindex, nullindex, nullindex } );
    data_edge_firstparent_triangle.resize    ( 2 * counter_edges + 3 * counter_triangles, nullindex                           );
    data_triangle_nextparents_of_edges.resize( 4 * counter_triangles                    , { nullindex, nullindex, nullindex } );
    
    data_triangle_vertices.resize               ( 4 * counter_triangles           , { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_triangle.resize     ( counter_vertices + counter_edges, nullindex                           );
    data_triangle_nextparents_of_vertices.resize( 4 * counter_triangles           , { nullindex, nullindex, nullindex } );
    
    data_edge_vertices.resize               ( 2 * counter_edges + 3 * counter_triangles, { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + counter_edges         , nullindex                );
    data_edge_nextparents_of_vertices.resize( 2 * counter_edges + 3 * counter_triangles, { nullindex, nullindex } );
    
    getcoordinates().addcoordinates( counter_edges );
    
    
    /*
     * update the flags of the 
     */
    
    flags_triangles.resize ( 4 * counter_triangles                    , SimplexFlag::SimplexFlagInvalid );
    flags_edges.resize     ( 2 * counter_edges + 3 * counter_triangles, SimplexFlag::SimplexFlagInvalid );
    flags_vertices.resize  ( counter_vertices + counter_edges         , SimplexFlag::SimplexFlagInvalid );
    
    // flags of newly created vertices
    for( int e = 0; e < counter_edges; e++ ) flags_vertices[ counter_vertices + e ] = flags_edges[e];
    
    // flags of bisected edges 
    for( int e = 0; e < counter_edges; e++ ) flags_edges[ counter_edges + e ] = flags_edges[e];
    
    // flags of ex nihilio new edges 
    for( int t = 0; t < counter_triangles; t++ ) {
        flags_edges[ 2 * counter_edges + 0 * counter_triangles + t ] = flags_triangles[t];
        flags_edges[ 2 * counter_edges + 1 * counter_triangles + t ] = flags_triangles[t];
        flags_edges[ 2 * counter_edges + 2 * counter_triangles + t ] = flags_triangles[t];
    }
    
    // flags of ex nihilio new triangles 
    for( int t = 0; t < counter_triangles; t++ ) {
        flags_triangles[ 1 * counter_triangles + t ] = flags_triangles[t];
        flags_triangles[ 2 * counter_triangles + t ] = flags_triangles[t];
        flags_triangles[ 3 * counter_triangles + t ] = flags_triangles[t];
    }
    
    
    
    
    
    /* create the new coordinates and fill them up */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      getcoordinates().loadvector( counter_vertices + e, get_edge_midpoint( e ) );
    }
    
    
    
    
    /**** Edge -- Vertex Relation ****/
    
    /* for each old vertex, set the new parent edge */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_edge[v];
      
      assert( p != nullindex && 0 <= p && p < counter_edges );
      
      int vi = data_edge_vertices[p][0] == v ? 0 : 1;
      
      assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
      assert( data_edge_vertices[p][vi] == v );
      
      data_vertex_firstparent_edge[v] = vi * counter_edges + p;
    }
    
    
    /* for each old edge, relocate the data of the old vertices' old parent edges */
    
    for( int e  = 0; e  < counter_edges;  e++ )
    for( int vi = 0; vi <=            1; vi++ )
    {
      int q = data_edge_nextparents_of_vertices[e][vi];
      
      int v = data_edge_vertices[e][vi];
      
      assert( v != nullindex && 0 <= v && v < counter_vertices );
      
      if( q == nullindex ) {
        
        data_edge_nextparents_of_vertices[ vi * counter_edges + e ][vi] = nullindex;
        
      } else {
        
        assert( 0 <= q && q < counter_edges );
        
        int vinp = data_edge_vertices[q][0] == v ? 0 : 1;
        
        assert( data_edge_vertices[q][0] == v || data_edge_vertices[q][1] == v );
        assert( data_edge_vertices[q][vinp] == v );
        
        data_edge_nextparents_of_vertices[ vi * counter_edges + e ][vi] = vinp * counter_edges + q;
      
      } 
      
    }
    
    
    /* for each new vertex, 
     * set the first and second parent edge from the old edge 
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      data_vertex_firstparent_edge[counter_vertices + e] = e;
      
      data_edge_nextparents_of_vertices[ 0 * counter_edges + e ][1] = counter_edges + e;
      data_edge_nextparents_of_vertices[ 1 * counter_edges + e ][0] = nullindex;
    }
    
    
    /* for each old edge, run over the adjacent triangles 
     * and add the corresponding new edges to the list of 
     * parent edges of new vertex.
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int t = data_edge_firstparent_triangle[e];
      
      while( t != nullindex ) {
        
        int ei   = nullindex;
        int e_1  = nullindex; 
        int e_2  = nullindex;
        int vi_1 = nullindex;
        int vi_2 = nullindex;
        
        // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
        
        if(        data_triangle_edges[t][0] == e ) {
          ei = 0; e_1 = 0; e_2 = 1; vi_1 = 0; vi_2 = 0; 
        } else if( data_triangle_edges[t][1] == e ) {
          ei = 1; e_1 = 0; e_2 = 2; vi_1 = 1; vi_2 = 0; 
        } else if( data_triangle_edges[t][2] == e ) {
          ei = 2; e_1 = 1; e_2 = 2; vi_1 = 1; vi_2 = 1; 
        } else {
          unreachable();
        }
        
        assert( ei  != nullindex && e_1 != nullindex && e_2 != nullindex );
        
        int old_first_parent = data_vertex_firstparent_edge[ counter_vertices + e ];
        
        data_vertex_firstparent_edge[ counter_vertices + e ]
          = 2 * counter_edges + e_1 * counter_triangles + t;
        
        data_edge_nextparents_of_vertices[ 2 * counter_edges + e_1 * counter_triangles + t ][ vi_1 ]
          = 2 * counter_edges + e_2 * counter_triangles + t;
        
        data_edge_nextparents_of_vertices[ 2 * counter_edges + e_2 * counter_triangles + t ][ vi_2 ]
          = old_first_parent;
        
        t = data_triangle_nextparents_of_edges[ t ][ ei ];
        
      }
      
    }
    
    /* for each edge created from an old edge, set the vertices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int vertex_back  = data_edge_vertices[e][0];
      int vertex_front = data_edge_vertices[e][1];
      
      data_edge_vertices[e                ][0] = vertex_back;
      data_edge_vertices[e                ][1] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][0] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][1] = vertex_front;
    }
    
    /* for each triangle, set the vertices of the new edge */
    
    for( int t = 0; t < counter_triangles; t++ )
    {
      data_edge_vertices[ 2 * counter_edges + 0 * counter_triangles + t ][0] = counter_vertices + data_triangle_edges[t][0];
      data_edge_vertices[ 2 * counter_edges + 0 * counter_triangles + t ][1] = counter_vertices + data_triangle_edges[t][1];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_triangles + t ][0] = counter_vertices + data_triangle_edges[t][0];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_triangles + t ][1] = counter_vertices + data_triangle_edges[t][2];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_triangles + t ][0] = counter_vertices + data_triangle_edges[t][1];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_triangles + t ][1] = counter_vertices + data_triangle_edges[t][2];
    }
    
    
    /**** Triangle -- Vertex Relation ****/
    
    /* for each old vertex, set the new parent triangle */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_triangle[v];
      
      assert( p != nullindex && 0 <= p && p < counter_triangles );
      
      int vi = data_triangle_vertices[p][0] == v ? 0 : data_triangle_vertices[p][1] == v ? 1 : 2;
      
      assert( data_triangle_vertices[p][0] == v || data_triangle_vertices[p][1] == v  || data_triangle_vertices[p][2] == v );
      assert( data_triangle_vertices[p][vi] == v );
      
      data_vertex_firstparent_triangle[v] = vi * counter_triangles + p;
    }
    
    
    /* for each old triangle, relocate the data of the old vertices' parent triangle */
    
    for( int t  = 0; t  < counter_triangles;  t++ )
    for( int vi = 0; vi <                 3; vi++ )
    {
      int q = data_triangle_nextparents_of_vertices[t][vi];
      
      int v = data_triangle_vertices[t][vi];
      
      if( q == nullindex ) {
        
        data_triangle_nextparents_of_vertices[ vi * counter_triangles + t ][vi] = nullindex;
        
      } else {
        
        int vinp = data_triangle_vertices[q][0] == v ? 0 : data_triangle_vertices[q][1] == v ? 1 : 2;
        
        assert( data_triangle_vertices[q][0] == v || data_triangle_vertices[q][1] == v || data_triangle_vertices[q][2] == v );
        assert( data_triangle_vertices[q][vinp] == v );
        
        data_triangle_nextparents_of_vertices[ vi * counter_triangles + t ][vi] = vinp * counter_triangles + q;
      
      } 
      
    }
    
    /* for each old edge, run over the adjacent triangles 
     * and add the corresponding new triangles to the list of 
     * parent triangles of new vertex.
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int t = data_edge_firstparent_triangle[e];
      
      while( t != nullindex ) {
        
        int ei   = nullindex;
        int t_1  = nullindex; 
        int t_2  = nullindex;
        int t_3  = nullindex;
        int vi_1 = nullindex;
        int vi_2 = nullindex;
        int vi_3 = nullindex;
        
        // [ 01 02 12 ] -> [ 01 02 ] [ 01 12 ] [ 02 12 ]
        // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
        if(        data_triangle_edges[t][0] == e ) {
          ei = 0; 
          t_1 = 0; t_2 = 3; t_3 = 1; vi_1 = 1; vi_2 = 0; vi_3 = 0;  
        } else if( data_triangle_edges[t][1] == e ) {
          ei = 1; 
          t_1 = 0; t_2 = 3; t_3 = 2; vi_1 = 2; vi_2 = 1; vi_3 = 0;  
        } else if( data_triangle_edges[t][2] == e ) {
          ei = 2; 
          t_1 = 1; t_2 = 3; t_3 = 2; vi_1 = 2; vi_2 = 2; vi_3 = 1; 
        } else {
          unreachable();
        }
        
        int old_first_parent = data_vertex_firstparent_triangle[ counter_vertices + e ];
        
        data_vertex_firstparent_triangle[ counter_vertices + e ]
          = t_1 * counter_triangles + t;
        
        if( t_1 != 0 ) assert( data_triangle_nextparents_of_vertices[ t_1 * counter_triangles + t ][ vi_1 ] == nullindex );
        data_triangle_nextparents_of_vertices[ t_1 * counter_triangles + t ][ vi_1 ]
          = t_2 * counter_triangles + t;
        
        if( t_2 != 0 ) assert( data_triangle_nextparents_of_vertices[ t_2 * counter_triangles + t ][ vi_2 ] == nullindex );
        data_triangle_nextparents_of_vertices[ t_2 * counter_triangles + t ][ vi_2 ]
          = t_3 * counter_triangles + t;
        
        if( t_3 != 0 ) assert( data_triangle_nextparents_of_vertices[ t_3 * counter_triangles + t ][ vi_3 ] == nullindex );
        data_triangle_nextparents_of_vertices[ t_3 * counter_triangles + t ][ vi_3 ]
          = old_first_parent;
        
        t = data_triangle_nextparents_of_edges[ t ][ ei ];
        
      }
      
    }
    
    
    
//     /* for each triangle, create the new vertices */
//     
//     for( int t = 0; t < counter_triangles; t++ )
//     {
//       int v00 = data_triangle_vertices[t][0];
//       int v11 = data_triangle_vertices[t][1];
//       int v22 = data_triangle_vertices[t][2];
//       
//       int v01 = counter_vertices + data_triangle_edges[t][0];
//       int v02 = counter_vertices + data_triangle_edges[t][1];
//       int v12 = counter_vertices + data_triangle_edges[t][2];
//       
//       // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
//         
//       data_triangle_vertices[ 0 * counter_triangles + t ][0] = v00;
//       data_triangle_vertices[ 0 * counter_triangles + t ][1] = v01;
//       data_triangle_vertices[ 0 * counter_triangles + t ][2] = v02;
//       
//       data_triangle_vertices[ 1 * counter_triangles + t ][0] = v01;
//       data_triangle_vertices[ 1 * counter_triangles + t ][1] = v11;
//       data_triangle_vertices[ 1 * counter_triangles + t ][2] = v12;
//       
//       data_triangle_vertices[ 2 * counter_triangles + t ][0] = v02;
//       data_triangle_vertices[ 2 * counter_triangles + t ][1] = v12;
//       data_triangle_vertices[ 2 * counter_triangles + t ][2] = v22;
//       
//       data_triangle_vertices[ 3 * counter_triangles + t ][0] = v01;
//       data_triangle_vertices[ 3 * counter_triangles + t ][1] = v02;
//       data_triangle_vertices[ 3 * counter_triangles + t ][2] = v12;
//       
//     }
    
    
    // TODO: Check code until here.
    /**** Triangle -- Edge Relation ****/
    
    /* for each old edge, set the new first parent triangle */
    /* for each new edge, set the new first parent triangle */
    // checked
    for( int e = 0; e < counter_edges; e++ )
    {
      int p = data_edge_firstparent_triangle[e];
      
      assert( p != nullindex );
      
      int ei        = nullindex;
      int nfp_back  = nullindex;
      int nfp_front = nullindex;
      
      if( data_triangle_edges[p][0] == e ){
        ei = 0; nfp_back = 0; nfp_front = 1;
      } else if( data_triangle_edges[p][1] == e ) {
        ei = 1; nfp_back = 0; nfp_front = 2;
      } else if( data_triangle_edges[p][2] == e ) {
        ei = 2; nfp_back = 1; nfp_front = 2;
      } else {
        unreachable();
      }
      
      assert( ei != nullindex );
      assert( data_triangle_edges[p][0] == e || data_triangle_edges[p][1] == e || data_triangle_edges[p][2] == e );
      assert( data_triangle_edges[p][ei] == e );
      
      data_edge_firstparent_triangle[ 0 * counter_edges + e ] = nfp_back  * counter_triangles + p;
      data_edge_firstparent_triangle[ 1 * counter_edges + e ] = nfp_front * counter_triangles + p;
      
    }
    
    
    /* for each triangle, relocate the data of the old edges' parent triangle */
    /* additionally, set the new parents */
    // checked
    for( int t  = 0; t  < counter_triangles;  t++ )
    for( int ei = 0; ei <                 3; ei++ )
    {
      int e = data_triangle_edges[t][ei];
      
      int t_back  = nullindex;
      int t_front = nullindex;
      int e_back  = nullindex;
      int e_front = nullindex;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
      
      if( ei == 0 ){
        t_back = 0; t_front = 1; e_back = 0; e_front = 0; 
      } else if( ei == 1 ) {
        t_back = 0; t_front = 2; e_back = 1; e_front = 1; 
      } else if( ei == 2 ) {
        t_back = 1; t_front = 2; e_back = 2; e_front = 2; 
      } else {
        unreachable();
      }
      
      
      int q = data_triangle_nextparents_of_edges[t][ei];
      
      if( q == nullindex ) {
        
        data_triangle_nextparents_of_edges[ t_back  * counter_triangles + t ][ e_back  ] = nullindex;
        data_triangle_nextparents_of_edges[ t_front * counter_triangles + t ][ e_front ] = nullindex;
        
      } else {
        
        int q_ei        = nullindex;
        int q_nfp_back  = nullindex;
        int q_nfp_front = nullindex;
        
        if(        data_triangle_edges[q][0] == e ){
          q_ei = 0; q_nfp_back = 0; q_nfp_front = 1;
        } else if( data_triangle_edges[q][1] == e ) {
          q_ei = 1; q_nfp_back = 0; q_nfp_front = 2;
        } else if( data_triangle_edges[q][2] == e ) {
          q_ei = 2; q_nfp_back = 1; q_nfp_front = 2;
        } else {
          unreachable();
        }
        
        assert( q_ei != nullindex );
        assert( data_triangle_edges[q][0] == e || data_triangle_edges[q][1] == e || data_triangle_edges[q][2] == e );
        assert( data_triangle_edges[q][q_ei] == e );
        
        data_triangle_nextparents_of_edges[ t_back  * counter_triangles + t ][ e_back  ] = q_nfp_back  * counter_triangles + q;
        data_triangle_nextparents_of_edges[ t_front * counter_triangles + t ][ e_front ] = q_nfp_front * counter_triangles + q;
        
      } 
      
    }
    
    /* for each triangle, run over the new edges and add firstparents and parents */
    // checked 
    for( int t = 0; t < counter_triangles; t++ )
    {
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
      // new edges [ 01 02 ], [ 01 12 ], [ 02 12 ]
      
      data_edge_firstparent_triangle[ 2 * counter_edges + 0 * counter_triangles + t ] = 3 * counter_triangles + t;
      data_edge_firstparent_triangle[ 2 * counter_edges + 1 * counter_triangles + t ] = 3 * counter_triangles + t;
      data_edge_firstparent_triangle[ 2 * counter_edges + 2 * counter_triangles + t ] = 3 * counter_triangles + t;
      
      data_triangle_nextparents_of_edges[ 3 * counter_triangles + t ][0] = 0 * counter_triangles + t;
      data_triangle_nextparents_of_edges[ 3 * counter_triangles + t ][1] = 1 * counter_triangles + t;
      data_triangle_nextparents_of_edges[ 3 * counter_triangles + t ][2] = 2 * counter_triangles + t;
      
      data_triangle_nextparents_of_edges[ 0 * counter_triangles + t ][2] = nullindex;
      data_triangle_nextparents_of_edges[ 1 * counter_triangles + t ][1] = nullindex;
      data_triangle_nextparents_of_edges[ 2 * counter_triangles + t ][0] = nullindex;
      
    }
    
    
    
    /* for each new triangle, create the new vertices */
    // checked
    for( int t = 0; t < counter_triangles; t++ )
    {
      int v00 = data_triangle_vertices[t][0];
      int v11 = data_triangle_vertices[t][1];
      int v22 = data_triangle_vertices[t][2];
      
      int v01 = counter_vertices + data_triangle_edges[t][0];
      int v02 = counter_vertices + data_triangle_edges[t][1];
      int v12 = counter_vertices + data_triangle_edges[t][2];
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_triangle_vertices[ 0 * counter_triangles + t ][0] = v00;
      data_triangle_vertices[ 0 * counter_triangles + t ][1] = v01;
      data_triangle_vertices[ 0 * counter_triangles + t ][2] = v02;
      
      data_triangle_vertices[ 1 * counter_triangles + t ][0] = v01;
      data_triangle_vertices[ 1 * counter_triangles + t ][1] = v11;
      data_triangle_vertices[ 1 * counter_triangles + t ][2] = v12;
      
      data_triangle_vertices[ 2 * counter_triangles + t ][0] = v02;
      data_triangle_vertices[ 2 * counter_triangles + t ][1] = v12;
      data_triangle_vertices[ 2 * counter_triangles + t ][2] = v22;
      
      data_triangle_vertices[ 3 * counter_triangles + t ][0] = v01;
      data_triangle_vertices[ 3 * counter_triangles + t ][1] = v02;
      data_triangle_vertices[ 3 * counter_triangles + t ][2] = v12;
      
    }
    
    /* for each new triangle, create the new edges */
    // checked
    for( int t = 0; t < counter_triangles; t++ )
    {
      int e00_01 = 0 * counter_edges + data_triangle_edges[t][0];
      int e00_02 = 0 * counter_edges + data_triangle_edges[t][1];
      int e01_11 = 1 * counter_edges + data_triangle_edges[t][0];
      int e02_22 = 1 * counter_edges + data_triangle_edges[t][1];
      int e11_12 = 0 * counter_edges + data_triangle_edges[t][2];
      int e12_22 = 1 * counter_edges + data_triangle_edges[t][2];
      
      int e01_02 = 2 * counter_edges + 0 * counter_triangles + t;
      int e01_12 = 2 * counter_edges + 1 * counter_triangles + t;
      int e02_12 = 2 * counter_edges + 2 * counter_triangles + t;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_triangle_edges[ 0 * counter_triangles + t ][0] = e00_01;
      data_triangle_edges[ 0 * counter_triangles + t ][1] = e00_02;
      data_triangle_edges[ 0 * counter_triangles + t ][2] = e01_02;
      
      data_triangle_edges[ 1 * counter_triangles + t ][0] = e01_11;
      data_triangle_edges[ 1 * counter_triangles + t ][1] = e01_12;
      data_triangle_edges[ 1 * counter_triangles + t ][2] = e11_12;
      
      data_triangle_edges[ 2 * counter_triangles + t ][0] = e02_12;
      data_triangle_edges[ 2 * counter_triangles + t ][1] = e02_22;
      data_triangle_edges[ 2 * counter_triangles + t ][2] = e12_22;
      
      data_triangle_edges[ 3 * counter_triangles + t ][0] = e01_02;
      data_triangle_edges[ 3 * counter_triangles + t ][1] = e01_12;
      data_triangle_edges[ 3 * counter_triangles + t ][2] = e02_12;
      
    }
    
    
    
    
    /* update the counters */
    
    counter_vertices  = counter_vertices + counter_edges;
    counter_edges     = 2 * counter_edges + 3 * counter_triangles;
    counter_triangles = 4 * counter_triangles;
    
    
    /* DONE */
    
    check();
}










void MeshSimplicial2D::midpoint_refinement( int t )
{
    // MeshSimplicial2D::check();
    
    assert( 0 <= t && t < counter_triangles );
    
    /* Allocate memory */
    
    data_triangle_edges.resize               ( counter_triangles + 2, { nullindex, nullindex, nullindex } );
    data_edge_firstparent_triangle.resize    ( counter_edges + 3,                             nullindex   );
    data_triangle_nextparents_of_edges.resize( counter_triangles + 2, { nullindex, nullindex, nullindex } );
    
    data_triangle_vertices.resize               ( counter_triangles + 2, { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_triangle.resize     ( counter_vertices + 1,                          nullindex   );
    data_triangle_nextparents_of_vertices.resize( counter_triangles + 2, { nullindex, nullindex, nullindex } );
    
    data_edge_vertices.resize               ( counter_edges + 3,   { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1,             nullindex   );
    data_edge_nextparents_of_vertices.resize( counter_edges + 3,   { nullindex, nullindex } );
    
    getcoordinates().addcoordinates( 1 );
    
    
    /* load the new coordinate */
    
    getcoordinates().loadvector( counter_vertices, get_triangle_midpoint( t ) );
    
    
    /* assemble the data and auxiliary variables */
    int t_e0 = data_triangle_edges[t][0];
    int t_e1 = data_triangle_edges[t][1];
    int t_e2 = data_triangle_edges[t][2];
    
    int t_v0 = data_triangle_vertices[t][0];
    int t_v1 = data_triangle_vertices[t][1];
    int t_v2 = data_triangle_vertices[t][2];
    
    int t_e_n0 = data_triangle_nextparents_of_edges[t][0];
    int t_e_n1 = data_triangle_nextparents_of_edges[t][1];
    int t_e_n2 = data_triangle_nextparents_of_edges[t][2];
    
    int t_v_n0 = data_triangle_nextparents_of_vertices[t][0];
    int t_v_n1 = data_triangle_nextparents_of_vertices[t][1];
    int t_v_n2 = data_triangle_nextparents_of_vertices[t][2];
    
    
    int t0 = t;
    int t1 = counter_triangles;
    int t2 = counter_triangles + 1;
    
    int e0 = counter_edges + 0; //  0 -> vn
    int e1 = counter_edges + 1; //  1 -> vn
    int e2 = counter_edges + 2; // vn -> 2
    
    int vn = counter_vertices;
    
    
    /* fill in: downward */
    data_triangle_vertices[ t0 ][0] = t_v0;
    data_triangle_vertices[ t0 ][1] = t_v1;
    data_triangle_vertices[ t0 ][2] = vn;
    
    data_triangle_vertices[ t1 ][0] = t_v0;
    data_triangle_vertices[ t1 ][1] = vn;
    data_triangle_vertices[ t1 ][2] = t_v2;
    
    data_triangle_vertices[ t2 ][0] = t_v1;
    data_triangle_vertices[ t2 ][1] = vn;
    data_triangle_vertices[ t2 ][2] = t_v2;
    
    data_triangle_edges[ t0 ][0] = t_e0;
    data_triangle_edges[ t0 ][1] = e0;
    data_triangle_edges[ t0 ][2] = e1;
    
    data_triangle_edges[ t1 ][0] = e0;
    data_triangle_edges[ t1 ][1] = t_e1;
    data_triangle_edges[ t1 ][2] = e2;
    
    data_triangle_edges[ t2 ][0] = e1;
    data_triangle_edges[ t2 ][1] = t_e2;
    data_triangle_edges[ t2 ][2] = e2;
    
    data_edge_vertices[ e0 ][0] = t_v0;
    data_edge_vertices[ e0 ][1] = vn;
    
    data_edge_vertices[ e1 ][0] = t_v1;
    data_edge_vertices[ e1 ][1] = vn;
    
    data_edge_vertices[ e2 ][0] = vn;
    data_edge_vertices[ e2 ][1] = t_v2;
    
    /* fill in: parentlist */
    data_triangle_nextparents_of_edges[ t0 ][0] = t_e_n0;
    data_triangle_nextparents_of_edges[ t0 ][1] = t1;
    data_triangle_nextparents_of_edges[ t0 ][2] = t2;
    
    data_triangle_nextparents_of_edges[ t1 ][0] = nullindex;
    data_triangle_nextparents_of_edges[ t1 ][1] = t_e_n1;
    data_triangle_nextparents_of_edges[ t1 ][2] = t2;
    
    data_triangle_nextparents_of_edges[ t2 ][0] = nullindex;
    data_triangle_nextparents_of_edges[ t2 ][1] = t_e_n2;
    data_triangle_nextparents_of_edges[ t2 ][2] = nullindex;
    
    data_triangle_nextparents_of_vertices[ t0 ][0] = t1;
    data_triangle_nextparents_of_vertices[ t0 ][1] = t2;
    data_triangle_nextparents_of_vertices[ t0 ][2] = t1;
    
    data_triangle_nextparents_of_vertices[ t1 ][0] = t_v_n0;
    data_triangle_nextparents_of_vertices[ t1 ][1] = t2;
    data_triangle_nextparents_of_vertices[ t1 ][2] = t2;
    
    data_triangle_nextparents_of_vertices[ t2 ][0] = t_v_n1;
    data_triangle_nextparents_of_vertices[ t2 ][1] = nullindex;
    data_triangle_nextparents_of_vertices[ t2 ][2] = t_v_n2;
    
    data_edge_nextparents_of_vertices[ e0 ][0] = data_vertex_firstparent_edge[ t_v0 ];
    data_edge_nextparents_of_vertices[ e0 ][1] = e1;
    data_vertex_firstparent_edge[ t_v0 ] = e0;
    
    data_edge_nextparents_of_vertices[ e1 ][0] = data_vertex_firstparent_edge[ t_v1 ];
    data_edge_nextparents_of_vertices[ e1 ][1] = e2;
    data_vertex_firstparent_edge[ t_v1 ] = e1;
    
    data_edge_nextparents_of_vertices[ e2 ][0] = nullindex;
    data_edge_nextparents_of_vertices[ e2 ][1] = data_vertex_firstparent_edge[ t_v2 ];
    data_vertex_firstparent_edge[ t_v2 ] = e2;
    
    /* edge t_e0: nothing needs to change */
    
    /* edge t_e1: relink */
    if( data_edge_firstparent_triangle[ t_e1 ] == t ) {
      data_edge_firstparent_triangle[ t_e1 ] = t1;
    } else {
      int current_triangle = data_edge_firstparent_triangle[ t_e1 ];
      while( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != t 
             &&
             data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != nullindex )
        current_triangle = data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ];
      assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] != nullindex );
      assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] == t );
      data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e1 ) ] = t1;
    }
        
    /* edge t_e2: relink */
    if( data_edge_firstparent_triangle[ t_e2 ] == t ) {
      data_edge_firstparent_triangle[ t_e2 ] = t2;
    } else {
      int current_triangle = data_edge_firstparent_triangle[ t_e2 ];
      while( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != t 
             &&
             data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex )
        current_triangle = data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ];
      assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] != nullindex );
      assert( data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] == t );
      data_triangle_nextparents_of_edges[ current_triangle ][ indexof_triangle_edge( current_triangle, t_e2 ) ] = t2;
    }
    
    /* vertex t_v0: nothing needs to change */
    
    /* vertex t_v1: nothing needs to change */
    
    /* vertex t_v2: relink */ 
    if( data_vertex_firstparent_triangle[ t_v2 ] == t ) {
      data_vertex_firstparent_triangle[ t_v2 ] = t1;
    } else {
      int current_triangle = data_vertex_firstparent_triangle[ t_v2 ];
      while( data_triangle_nextparents_of_vertices[ current_triangle ][ indexof_triangle_vertex( current_triangle, t_v2 ) ] != t 
             &&
             data_triangle_nextparents_of_vertices[ current_triangle ][ indexof_triangle_vertex( current_triangle, t_v2 ) ] != nullindex )
        current_triangle = data_triangle_nextparents_of_vertices[ current_triangle ][ indexof_triangle_vertex( current_triangle, t_v2 ) ];
      assert( data_triangle_nextparents_of_vertices[ current_triangle ][ indexof_triangle_vertex( current_triangle, t_v2 ) ] != nullindex );
      assert( data_triangle_nextparents_of_vertices[ current_triangle ][ indexof_triangle_vertex( current_triangle, t_v2 ) ] == t );
      data_triangle_nextparents_of_vertices[ current_triangle ][ indexof_triangle_vertex( current_triangle, t_v2 ) ] = t1;
    }
    
    /* edge e0: link */
    data_edge_firstparent_triangle[ e0 ] = t0;
    
    /* edge e1: link */
    data_edge_firstparent_triangle[ e1 ] = t0;
    
    /* edge e2: link */
    data_edge_firstparent_triangle[ e2 ] = t1;
    
    /* vertex vn: link */
    data_vertex_firstparent_triangle[ counter_vertices ] = t0;
    data_vertex_firstparent_edge[ counter_vertices ] = e0;
    
    
    
    /*
     * update the flags of the 
     */
    
    flags_triangles.resize ( counter_triangles + 2, SimplexFlag::SimplexFlagInvalid );
    flags_edges.resize     ( counter_edges + 3    , SimplexFlag::SimplexFlagInvalid );
    flags_vertices.resize  ( counter_vertices + 1 , SimplexFlag::SimplexFlagInvalid );
    
    flags_vertices[ counter_vertices ] = flags_triangles[t];
    
    flags_edges[ counter_edges + 0 ] = flags_triangles[t];
    flags_edges[ counter_edges + 1 ] = flags_triangles[t];
    flags_edges[ counter_edges + 2 ] = flags_triangles[t];

    flags_triangles[ counter_triangles + 0 ] = flags_triangles[t];
    flags_triangles[ counter_triangles + 1 ] = flags_triangles[t];
    
        
    
    
    
    /* Counters */
    counter_triangles = counter_triangles + 2;
    counter_edges     = counter_edges     + 3;
    counter_vertices  = counter_vertices  + 1;
    
    /* DONE */
    
    // MeshSimplicial2D::check();
}



void MeshSimplicial2D::midpoint_refinement_global()
{
    check();
    
    int N = counter_triangles;
    
    for( int t = 0; t < N; t++ ) {
      LOG << t << nl;
      midpoint_refinement( t );
      
    }
    
    check();
}
















FloatVector MeshSimplicial2D::get_triangle_midpoint( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_triangle_vertices(t)[0], d ) 
                 + getcoordinates().getdata( get_triangle_vertices(t)[1], d ) 
                 + getcoordinates().getdata( get_triangle_vertices(t)[2], d )
                ) / 3.;
    return mid;
}

FloatVector MeshSimplicial2D::get_edge_midpoint    ( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_edge_vertices(e)[0], d ) 
                 + getcoordinates().getdata( get_edge_vertices(e)[1], d )
                ) / 2.;
    return mid;
}

Float MeshSimplicial2D::get_edge_length( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    Float length = 0.;
    for( int d = 0; d < getouterdimension(); d++ )
      length += power_numerical( getcoordinates().getdata( get_edge_vertices(e)[0], d ) - getcoordinates().getdata( get_edge_vertices(e)[1], d ), 2. );
    return std::sqrt( length );
}
        

int MeshSimplicial2D::get_oldest_edge( int t ) const
{
    assert( 0 <= t && t < counter_triangles );
    int e0 = get_triangle_edge( t, 0 );
    int e1 = get_triangle_edge( t, 1 );
    int e2 = get_triangle_edge( t, 2 );
    if( e0 < e1 and e0 < e2 ) return e0;
    if( e1 < e0 and e1 < e2 ) return e1;
    if( e2 < e0 and e2 < e1 ) return e2;
    unreachable();
}










void MeshSimplicial2D::merge( const MeshSimplicial2D& mesh )
{
    check();
    mesh.check();
    
    // Edges -> vertices 
    
    {
        auto& mine         =      data_edge_vertices;
        const auto& theirs = mesh.data_edge_vertices;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ][0] = theirs[i][0] + counter_vertices;
            mine[ old_size + i ][1] = theirs[i][1] + counter_vertices;
        }
    }
    
    {
        auto& mine         =      data_vertex_firstparent_edge;
        const auto& theirs = mesh.data_vertex_firstparent_edge;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ] = ( theirs[i] == nullindex ) ? ( nullindex ) : theirs[i] + counter_edges;
        }
    }
    
    {
        auto& mine         =      data_edge_nextparents_of_vertices;
        const auto& theirs = mesh.data_edge_nextparents_of_vertices;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ][0] = ( theirs[i][0] == nullindex ) ? ( nullindex ) : theirs[i][0] + counter_edges;
            mine[ old_size + i ][1] = ( theirs[i][1] == nullindex ) ? ( nullindex ) : theirs[i][1] + counter_edges;
        }
    }
    
    
    
    // Triangles -> vertices 
    
    {
        auto& mine         =      data_triangle_vertices;
        const auto& theirs = mesh.data_triangle_vertices;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ][0] = theirs[i][0] + counter_vertices;
            mine[ old_size + i ][1] = theirs[i][1] + counter_vertices;
            mine[ old_size + i ][2] = theirs[i][2] + counter_vertices;
        }
    }
    
    {
        auto& mine         =      data_vertex_firstparent_triangle;
        const auto& theirs = mesh.data_vertex_firstparent_triangle;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ] = theirs[i] + counter_triangles;
        }
    }
    
    {
        auto& mine         =      data_triangle_nextparents_of_vertices;
        const auto& theirs = mesh.data_triangle_nextparents_of_vertices;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ][0] = ( theirs[i][0] == nullindex ) ? ( nullindex ) : theirs[i][0] + counter_triangles;
            mine[ old_size + i ][1] = ( theirs[i][1] == nullindex ) ? ( nullindex ) : theirs[i][1] + counter_triangles;
            mine[ old_size + i ][2] = ( theirs[i][2] == nullindex ) ? ( nullindex ) : theirs[i][2] + counter_triangles;
        }
    }
    
    
    
    // Triangles -> edges 
    
    {
        auto& mine         =      data_triangle_edges;
        const auto& theirs = mesh.data_triangle_edges;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ][0] = theirs[i][0] + counter_edges;
            mine[ old_size + i ][1] = theirs[i][1] + counter_edges;
            mine[ old_size + i ][2] = theirs[i][2] + counter_edges;
        }
    }
    
    {
        auto& mine         =      data_edge_firstparent_triangle;
        const auto& theirs = mesh.data_edge_firstparent_triangle;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ] = theirs[i] + counter_triangles;
        }
    }
    
    {
        auto& mine         =      data_triangle_nextparents_of_edges;
        const auto& theirs = mesh.data_triangle_nextparents_of_edges;
        auto old_size = mine.size();
        
        mine.resize( mine.size() + theirs.size() );
        for( int i = 0; i < theirs.size(); i++ ) {
            mine[ old_size + i ][0] = ( theirs[i][0] == nullindex ) ? ( nullindex ) : theirs[i][0] + counter_triangles;
            mine[ old_size + i ][1] = ( theirs[i][1] == nullindex ) ? ( nullindex ) : theirs[i][1] + counter_triangles;
            mine[ old_size + i ][2] = ( theirs[i][2] == nullindex ) ? ( nullindex ) : theirs[i][2] + counter_triangles;
        }
    }
    
    
    
    // Simplex flags 
    
    {
        auto& mine         =      flags_vertices;
        const auto& theirs = mesh.flags_vertices;
        mine.insert( mine.end(), theirs.begin(), theirs.end() );
    }
    
    {
        auto& mine         =      flags_edges;
        const auto& theirs = mesh.flags_edges;
        mine.insert( mine.end(), theirs.begin(), theirs.end() );
    }
    
    {
        auto& mine         =      flags_triangles;
        const auto& theirs = mesh.flags_triangles;
        mine.insert( mine.end(), theirs.begin(), theirs.end() );
    }
    
    
    // counters 
    
    counter_vertices  += mesh.counter_vertices;
    counter_edges     += mesh.counter_edges;
    counter_triangles += mesh.counter_triangles;
    
    getcoordinates().append( mesh.getcoordinates() );
    
    check();
}




















std::string MeshSimplicial2D::outputTikZ() const
{
    std::ostringstream os;
    
    const auto& coords = getcoordinates();
    
    for( int n = 0; n < coords.getnumber(); n++ )
    {
        os << "\\coordinate ("
           << n
           << ")  at ( ";
        for( int d = 0; d < coords.getdimension(); d++ )
        {
           os << coords.getdata(n,d) 
              << ( d == coords.getdimension()-1 ? "" : ", ");
        }
        
        os << ");" << nl;
        
    }
    
    for( int t = 0; t < count_simplices(2); t++ )
    {
        os << "\\draw ";
        os << "(" << get_subsimplex(2,0,t,0) << ") -- ";
        os << "(" << get_subsimplex(2,0,t,1) << ") -- ";
        os << "(" << get_subsimplex(2,0,t,2) << ") -- ";
        os << "cycle;" << nl;
    }
    
    return os.str();
}
        


inline std::string rgb_to_string( unsigned char r, unsigned char g, unsigned char b )
{
    char result[8] = {0,0,0,0,0,0,0,0};
    snprintf(result, countof(result), "#%02x%02x%02x", r, g, b);
    return std::string( result );
}

inline unsigned char float_to_color( Float f )
{
  if( f >= 255.) f = 255.;
  if( f <= 0.  ) f = 0.;
  return (unsigned char)f;
}

UNUSED inline int leading_digits( double num )
{
    // If the number is negative, take its absolute value
    // Calculate the order of magnitude of the number
    // If the number is greater than 1, return the integer part of the order, else, return 0
    
    if( num < 0 ) num = -num;
    int order = static_cast<int>( ceil( log10(num) ) );
    if( order >= 0 ) return order + 1; else return 1;
}

inline std::string render_number( double num, unsigned int tail = 8 )
{
    assert( std::isfinite(num) );
    
    // int lead = leading_digits( num );
    // const int str_number_of_chars = 1+lead+1+tail+1;
    // char str[str_number_of_chars];
    // snprintf( str, str_number_of_chars, "% *.*f", 1+lead+1+tail, tail, num);

    const int str_number_of_printed_chars = ( 1 + 1 + 1 ) + tail + 1 + (4);
    assert( str_number_of_printed_chars >= 0 );
    char str[str_number_of_printed_chars+1];
    snprintf( str, str_number_of_printed_chars, "%.*e", tail, num );
    
    return std::string(str);
}

std::string MeshSimplicial2D::outputSVG( 
    Float stroke_width,
    const std::string& fill,
    const std::string& stroke,
    const FloatVector* triangle_red,
    const FloatVector* triangle_green,
    const FloatVector* triangle_blue
) const {
    std::ostringstream os;

    auto coords = getcoordinates();

    coords.shift( FloatVector{ 
        -getcoordinates().getmin(0), 
        -getcoordinates().getmin(1) 
    } );

    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.200000\" width=\"100%\" height=\"100%\" "
       << "viewBox=\""
       << coords.getmin(0) << space << coords.getmin(1) << space 
       << coords.getmax(0) << space << coords.getmax(1) << "\""
       << " shape-rendering=\"crispEdges\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
       << nl;

    for( int t = 0; t < this->count_triangles(); t++ ) {
      
      int v0 = this->get_triangle_vertex(t,0);
      int v1 = this->get_triangle_vertex(t,1);
      int v2 = this->get_triangle_vertex(t,2);
      
      Float x0 = coords.getdata(v0,0);
      Float y0 = coords.getdata(v0,1);
      Float x1 = coords.getdata(v1,0);
      Float y1 = coords.getdata(v1,1);
      Float x2 = coords.getdata(v2,0);
      Float y2 = coords.getdata(v2,1);
      
      os << "<polygon points=\"" 
         << render_number(x0) << "," << render_number(y0) << " " 
         << render_number(x1) << "," << render_number(y1) << " " 
         << render_number(x2) << "," << render_number(y2) << "\"";
      
      if( fill == "array" ){
      
        assert( triangle_red ); assert( triangle_green ); assert( triangle_blue );
      
        double channel_red   = triangle_red->at(t);
        double channel_green = triangle_green->at(t);
        double channel_blue  = triangle_blue->at(t);

        // LOG << channel_blue << space << channel_green << space << channel_red << nl;

        assert( 0. <= channel_red   and channel_red   <= 255. );
        assert( 0. <= channel_green and channel_green <= 255. );
        assert( 0. <= channel_blue  and channel_blue  <= 255. );
      
        std::string colorstring = rgb_to_string( 
            (unsigned char)channel_red, 
            (unsigned char)channel_green, 
            (unsigned char)channel_blue );

        os << " fill=\"" << colorstring << "\"";
      } else {
        os << " fill=\"" << "red"       << "\"";
      }
      os << " stroke=\"" << stroke << "\"";
      os << " stroke-width=\"" << render_number(stroke_width) << "\"";
      os << "></polygon>";
      os << nl;
        
    }
    
    os << "</svg>";
    
    return os.str();
}

inline void computeProjection( Float x0, Float y0, Float x1, Float y1, Float x2, Float y2, Float& px, Float& py ) 
{
    Float dx = x2 - x1; Float dy = y2 - y1;
    Float norm_sq = dx * dx + dy * dy;
    Float lambda = ( (x0 - x1) * dx + (y0 - y1) * dy ) / norm_sq;
    px = x1 + lambda * dx; 
    py = y1 + lambda * dy;
}

std::string MeshSimplicial2D::outputLinearSVG( 
    const FloatVector& triangle_red,
    const FloatVector& triangle_green,
    const FloatVector& triangle_blue,
    Float stroke_width,
    const std::string& fill,
    const std::string& stroke
) const {
    std::ostringstream os;

    // 1. copy coordinates and shift them 
    auto coords = getcoordinates();

    coords.shift( FloatVector{ 
        -getcoordinates().getmin(0), 
        -getcoordinates().getmin(1) 
    } );

    // 2. preamble 
    os << "<svg xmlns=\"http://www.w3.org/2000/svg\" version=\"1.200000\" width=\"100%\" height=\"100%\" "
       << "viewBox=\""
       << coords.getmin(0) << space << coords.getmin(1) << space 
       << coords.getmax(0) << space << coords.getmax(1) << "\""
       << " shape-rendering=\"crispEdges\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">"
       << nl;

    // 2. print the defs of the linear gradients 
    os << "<defs>" << nl;
    for( int t = 0; t < this->count_triangles(); t++ ) {
      
      int v0 = this->get_triangle_vertex(t,0);
      int v1 = this->get_triangle_vertex(t,1);
      int v2 = this->get_triangle_vertex(t,2);
      
      Float x0 = coords.getdata(v0,0); Float y0 = coords.getdata(v0,1);
      Float x1 = coords.getdata(v1,0); Float y1 = coords.getdata(v1,1);
      Float x2 = coords.getdata(v2,0); Float y2 = coords.getdata(v2,1);

      Float cx0, cx1, cx2, cy0, cy1, cy2;
      computeProjection( x0, y0, x1, y1, x2, y2, cx0, cy0 );
      computeProjection( x1, y1, x0, y0, x2, y2, cx1, cy1 );
      computeProjection( x2, y2, x1, y1, x0, y0, cx2, cy2 );

      // Float cx0 = ( x1 + x2 ) / 2.;
      // Float cy0 = ( y1 + y2 ) / 2.;
      // Float cx1 = ( x0 + x2 ) / 2.;
      // Float cy1 = ( y0 + y2 ) / 2.;
      // Float cx2 = ( x0 + x1 ) / 2.;
      // Float cy2 = ( y0 + y1 ) / 2.;

      //std::string rgb_to_string( unsigned char r, unsigned char g, unsigned char b )

      unsigned char r0 = float_to_color( triangle_red  [3*t+0] );
      unsigned char g0 = float_to_color( triangle_green[3*t+0] );
      unsigned char b0 = float_to_color( triangle_blue [3*t+0] );
      unsigned char r1 = float_to_color( triangle_red  [3*t+1] );
      unsigned char g1 = float_to_color( triangle_green[3*t+1] );
      unsigned char b1 = float_to_color( triangle_blue [3*t+1] );
      unsigned char r2 = float_to_color( triangle_red  [3*t+2] );
      unsigned char g2 = float_to_color( triangle_green[3*t+2] );
      unsigned char b2 = float_to_color( triangle_blue [3*t+2] );
      
      std::string color0 = rgb_to_string( r0, g0, b0 );
      std::string color1 = rgb_to_string( r1, g1, b1 );
      std::string color2 = rgb_to_string( r2, g2, b2 );
      
      std::ostringstream grad0;
      std::ostringstream grad1;
      std::ostringstream grad2;

      grad0 << "<linearGradient gradientUnits=\"userSpaceOnUse\" id=\"myshade0_" << t << "\" x1=\"" << render_number(x0) << "\" y1=\"" << render_number(y0) << "\" x2=\"" << render_number(cx0) << "\" y2=\"" << render_number(cy0) << "\" >"
            << "<stop offset=\"0%\" stop-color=\"" << color0 << "\" stop-opacity=\"100%\" />"
            << "<stop offset=\"100%\" stop-color=\"" << color0 << "\" stop-opacity=\"100%\" />"
            << "</linearGradient>";
      grad1 << "<linearGradient gradientUnits=\"userSpaceOnUse\" id=\"myshade1_" << t << "\" x1=\"" << render_number(x1) << "\" y1=\"" << render_number(y1) << "\" x2=\"" << render_number(cx1) << "\" y2=\"" << render_number(cy1) << "\" >"
            << "<stop offset=\"0%\" stop-color=\"" << color1 << "\" stop-opacity=\"100%\" />"
            << "<stop offset=\"100%\" stop-color=\"" << color1 << "\" stop-opacity=\"0%\" />"
            << "</linearGradient>";
      grad2 << "<linearGradient gradientUnits=\"userSpaceOnUse\" id=\"myshade2_" << t << "\" x1=\"" << render_number(x2) << "\" y1=\"" << render_number(y2) << "\" x2=\"" << render_number(cx2) << "\" y2=\"" << render_number(cy2) << "\" >"
            << "<stop offset=\"0%\" stop-color=\"" << color2 << "\" stop-opacity=\"100%\" />"
            << "<stop offset=\"100%\" stop-color=\"" << color2 << "\" stop-opacity=\"0%\" />"
            << "</linearGradient>";
      
      os << grad0.str() << nl;
      os << grad1.str() << nl;
      os << grad2.str() << nl;
    }
    os << "</defs>" << nl;
    
    
    // 3. print the triangles 
    for( int t = 0; t < this->count_triangles(); t++ ) {
      
      int v0 = this->get_triangle_vertex(t,0);
      int v1 = this->get_triangle_vertex(t,1);
      int v2 = this->get_triangle_vertex(t,2);
      
      Float x0 = coords.getdata(v0,0); Float y0 = coords.getdata(v0,1);
      Float x1 = coords.getdata(v1,0); Float y1 = coords.getdata(v1,1);
      Float x2 = coords.getdata(v2,0); Float y2 = coords.getdata(v2,1);
      
      std::ostringstream tri_start; tri_start << "<polygon points=\"" << render_number(x0) << "," << render_number(y0) << " " << render_number(x1) << "," << render_number(y1) << " " << render_number(x2) << "," << render_number(y2) << "\"";
      std::ostringstream tri_fill0; tri_fill0 << " fill=\"url(#myshade0_" << t << ")\" ";
      std::ostringstream tri_fill1; tri_fill1 << " fill=\"url(#myshade1_" << t << ")\" ";
      std::ostringstream tri_fill2; tri_fill2 << " fill=\"url(#myshade2_" << t << ")\" ";
      std::ostringstream tri_end  ; tri_end   << " stroke=\"" << stroke << "\" stroke-width=\"" << render_number(stroke_width) << "\"></polygon>";
      
      os << tri_start.str() << tri_fill0.str() << tri_end.str() << nl;
      os << tri_start.str() << tri_fill1.str() << tri_end.str() << nl;
      os << tri_start.str() << tri_fill2.str() << tri_end.str() << nl;
        
    }
    
    // 4. close the SVG 
    os << "</svg>";
    
    return os.str();
}






std::size_t MeshSimplicial2D::memorysize() const
{
    std::size_t ret = 0;

    ret += getcoordinates().memorysize();

    ret += sizeof(getinnerdimension());
    ret += sizeof(getouterdimension());

    ret += sizeof( counter_triangles );
    ret += sizeof( counter_edges     );
    ret += sizeof( counter_vertices  );

    { const auto& D = data_triangle_edges;                   ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = data_edge_firstparent_triangle;        ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = data_triangle_nextparents_of_edges;    ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    
    { const auto& D = data_triangle_vertices;                ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = data_vertex_firstparent_triangle;      ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = data_triangle_nextparents_of_vertices; ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    
    { const auto& D = data_edge_vertices;                    ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = data_vertex_firstparent_edge;          ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = data_edge_nextparents_of_vertices;     ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    
    { const auto& D = flags_triangles; ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = flags_edges;     ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = flags_vertices;  ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };

    return ret;
}




