
#include <algorithm>
#include <sstream>
#include <map>
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

#include "mesh.simplicial1D.hpp"




MeshSimplicial1D::MeshSimplicial1D( int outerdim )
:
    Mesh( 1, outerdim ),
    
    counter_edges(0),
    counter_vertices(0),
    
    data_edge_vertices(0),
    data_vertex_firstparent_edge(0),
    data_edge_nextparents_of_vertices(0),
    
    flags_edges   ( 0, SimplexFlag::SimplexFlagNull ),
    flags_vertices( 0, SimplexFlag::SimplexFlagNull )
    
{
    MeshSimplicial1D::check();
}


MeshSimplicial1D::MeshSimplicial1D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,2>>& edge_vertices
)
:
    Mesh( 1, outerdim ),
    
    counter_edges( edge_vertices.size() ),
    counter_vertices(0),
    
    data_edge_vertices( edge_vertices ),
    data_vertex_firstparent_edge( 0 ),
    data_edge_nextparents_of_vertices( counter_edges, { nullindex, nullindex } ),

    flags_edges   ( counter_edges, SimplexFlag::SimplexFlagNull ),
    flags_vertices( 0, SimplexFlag::SimplexFlagNull )
    
{
    
    getcoordinates() = coords;
    
    /* 1. Count edges, transfer data */ 
    /* DONE */
    
    
    /* 2. Count vertices, allocate memory */
    counter_vertices = 0;
    for( const auto& duple : data_edge_vertices )
    for( const int& vertex : duple )
      counter_vertices = counter_vertices < vertex ? vertex : counter_vertices; 
    counter_vertices += 1;
    
    data_vertex_firstparent_edge.resize( counter_vertices, nullindex );
    
    /* 3. For each vertex, set the first parent and the neighboring parents */
    
    for( int e =  0; e  <  counter_edges; e++  )
    for( int vi = 0; vi <=             1; vi++ )
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
    
    
    /* Coda: set the flags to the Null flag */
    flags_vertices.resize( counter_vertices, SimplexFlag::SimplexFlagInvalid );
    
    for( int e =  0; e  <  counter_edges;    e++  ) flags_edges.at( e )    = SimplexFlag::SimplexFlagNull;
    for( int v =  0; v  <  counter_vertices; v++  ) flags_vertices.at( v ) = SimplexFlag::SimplexFlagNull;
    
    
    MeshSimplicial1D::check();
}


MeshSimplicial1D::MeshSimplicial1D( 
    int outerdim,
    const Coordinates& coords,
    const std::vector<std::array<int,2>>& edge_vertices,
    const std::vector<std::array<int,2>>& edge_nextparent_of_vertices,
    const std::vector<int              >& vertex_firstparent_edge
)
:
    Mesh( 1, outerdim ),
    
    counter_edges( edge_vertices.size() ),
    counter_vertices( vertex_firstparent_edge.size() ),
    
    data_edge_vertices( edge_vertices ),
    data_vertex_firstparent_edge( vertex_firstparent_edge ),
    data_edge_nextparents_of_vertices( edge_nextparent_of_vertices ),

    flags_edges   ( counter_edges,    SimplexFlag::SimplexFlagInvalid ),
    flags_vertices( counter_vertices, SimplexFlag::SimplexFlagInvalid )
    
{
    
    getcoordinates() = coords;
    
    /* set the flags to the Null flag */    
    for( int e =  0; e  <  counter_edges;    e++  ) flags_edges.at( e )    = SimplexFlag::SimplexFlagNull;
    for( int v =  0; v  <  counter_vertices; v++  ) flags_vertices.at( v ) = SimplexFlag::SimplexFlagNull;

    MeshSimplicial1D::check();
}


MeshSimplicial1D::~MeshSimplicial1D()
{
    MeshSimplicial1D::check();
}



bool MeshSimplicial1D::is_equal_to( const MeshSimplicial1D& mesh ) const 
{
  return counter_edges == mesh.counter_edges
         &&
         counter_vertices == mesh.counter_vertices
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
         flags_edges == mesh.flags_edges
         &&
         flags_vertices == mesh.flags_vertices
         &&
         true;
}


void MeshSimplicial1D::check() const
{
    
    #ifdef NDEBUG
    return;
    #else 
    
    /* 1. Check the array sizes */
    
    assert( counter_edges == data_edge_vertices.size() );
    assert( counter_edges == data_edge_nextparents_of_vertices.size() );
    assert( counter_vertices == data_vertex_firstparent_edge.size() );
    
    assert( count_vertices() == getcoordinates().getnumber() );
    
    assert( counter_edges    == flags_edges.size()    );
    assert( counter_vertices == flags_vertices.size() );
    
    
    /* 
     * each edge: each vertex is a valid index
     * each edge: each vertex is unique 
     * each edge: the next parents make sense 
     */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        assert( data_edge_vertices[e][0] != nullindex );
        assert( data_edge_vertices[e][1] != nullindex );
        assert( 0 <= data_edge_vertices[e][0] && data_edge_vertices[e][0] < counter_vertices );
        assert( 0 <= data_edge_vertices[e][1] && data_edge_vertices[e][1] < counter_vertices );
        assert( data_edge_vertices[e][0] != data_edge_vertices[e][1] );
        
        if( data_edge_nextparents_of_vertices[e][0] != nullindex || data_edge_nextparents_of_vertices[e][1] != nullindex )
          assert( data_edge_nextparents_of_vertices[e][0] != data_edge_nextparents_of_vertices[e][1] );
        
        for( int vi = 0; vi < 2; vi++ ) {
          
          if( data_edge_nextparents_of_vertices[e][vi] != nullindex ) {
            
            assert( 0 <= data_edge_nextparents_of_vertices[e][vi] && data_edge_nextparents_of_vertices[e][vi] < counter_edges );
            
            assert( data_edge_vertices[ data_edge_nextparents_of_vertices[e][vi] ][0] == data_edge_vertices[e][vi] 
                    ||
                    data_edge_vertices[ data_edge_nextparents_of_vertices[e][vi] ][1] == data_edge_vertices[e][vi] 
                  );
                  
          }
          
        }
          
    }
    
    /* 
     * check that all edges are unique, even up to permutation
     */
    
    for( int e1 = 0; e1 < counter_edges; e1++ )
    for( int e2 = 0; e2 < counter_edges; e2++ )
    {
        if( e1 == e2 ) continue;
        assert( data_edge_vertices[e1][0] != data_edge_vertices[e2][0] || data_edge_vertices[e1][1] != data_edge_vertices[e2][1] );
        assert( data_edge_vertices[e1][0] != data_edge_vertices[e2][1] || data_edge_vertices[e1][1] != data_edge_vertices[e2][0] );
    }
    
    /* 
     * each vertex_firstparent: first parent is non-null 
     */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
        int p = data_vertex_firstparent_edge[v];
        
        assert( p != nullindex );
        assert( 0 <= p && p < counter_edges );
        
        assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
    }
    
    /* 
     * check that each is listed as parent somewhere 
     */
    
    for( int e  = 0; e  <  counter_edges; e++ )
    for( int vi = 0; vi <=             1; vi++ )
    {
      
      int v = data_edge_vertices[e][vi];
      
      assert( 0 <= v && v < counter_vertices );
      
      int p = data_vertex_firstparent_edge[v];
      
      assert( p != nullindex && 0 <= p && p < counter_edges );
      
      while( p != e && p != nullindex )
        p = data_edge_nextparents_of_vertices[p][ ( data_edge_vertices[p][0] == v ) ? 0 : 1 ];
      
      assert( p == e );
      
    }
    
    
    /*
     * check that all the flags are valid
     */
    
    for( int e  = 0; e  <  counter_edges; e++ )
        assert( flags_edges[e] != SimplexFlag::SimplexFlagInvalid );

    for( int v  = 0; v  <  counter_vertices; v++ )
        assert( flags_vertices[v] != SimplexFlag::SimplexFlagInvalid );
    
    Mesh::check();
    
    #endif
}






std::string MeshSimplicial1D::text() const
{
    std::ostringstream os;

    os << "Triangulation of 1D Manifold!" << nl;
    
    os << counter_edges << space << counter_vertices << nl;
    
    os << "Edge vertices" << nl;
    
    for( const auto& duple : data_edge_vertices )
      os << duple[0] << space << duple[1] << nl;
    
    os << "Vertex first parents" << nl;
    
    for( int fp : data_vertex_firstparent_edge )
      os << fp << nl;
    
    os << "Edge next parents " << nl;
    
    for( const auto& duple : data_edge_nextparents_of_vertices )
      os << duple[0] << space << duple[1] << nl;
    
    os << "Finished printing" << nl;
    
    return os.str();
}






bool MeshSimplicial1D::has_dimension_counted( int dim ) const
{
    assert( 0 <= dim && dim <= 1 );
    return true;
}

int MeshSimplicial1D::count_simplices( int dim ) const
{
  if( dim == 0 )
    return count_vertices();
  else if( dim == 1 )
    return count_edges();
  else
    unreachable();
}

bool MeshSimplicial1D::has_subsimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= 1 );
    return true;
}

IndexMap MeshSimplicial1D::getsubsimplices( int sup, int sub, int cell ) const
{
  assert( 0 <= sub && sub <= sup && sup <= 1 );
  assert( 0 <= cell );
  if( sup == 0 ) assert( cell <= count_vertices() );
  if( sup == 1 ) assert( cell <= count_edges()    );
  
  
  if( sup == 0 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_vertices()-1), { cell } );
    
  } else if( sup == 1 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    auto temp = get_edge_vertices(cell); 
    return IndexMap( IndexRange(0,1), IndexRange( 0, count_vertices()-1 ), std::vector<int>( temp.begin(), temp.end() ) );
    
  } else if( sup == 1 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    return IndexMap( IndexRange(0,0), IndexRange(0,count_edges()-1), { cell } );
    
  } else {
      
    unreachable(); 
    
  }
  
}

bool MeshSimplicial1D::has_supersimplices_listed( int sup, int sub ) const
{
    assert( 0 <= sub && sub <= sup && sup <= 1 );
    return true;
}

const std::vector<int> MeshSimplicial1D::getsupersimplices( int sup, int sub, int cell ) const
{
  if( sup == 0 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    return { cell };
    
  } else if( sup == 1 && sub == 0 ) {
    
    assert( 0 <= cell && cell < count_vertices() );
    auto temp = get_edge_parents_of_vertex( cell ); 
    return temp; //return std::vector<int>( temp.begin(), temp.end() );
    
  } else if( sup == 1 && sub == 1 ) {
    
    assert( 0 <= cell && cell < count_edges() );
    return { cell };
    
  } else {
    
    unreachable(); 
    
  }
  
}



SimplexFlag MeshSimplicial1D::get_flag( int dim, int cell ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    if( dim == 0 ) {
        assert( 0 <= cell && cell < count_vertices() );
        return flags_vertices[cell];
    } else if( dim == 1 ) {
        assert( 0 <= cell && cell < count_edges() );
        return flags_edges[cell];
    } else {
        unreachable();
    }
}
        
void MeshSimplicial1D::set_flag( int dim, int cell, SimplexFlag flag )
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    if( dim == 0 ) {
        assert( 0 <= cell && cell < count_vertices() );
        flags_vertices[cell] = flag;
    } else if( dim == 1 ) {
        assert( 0 <= cell && cell < count_edges() );
        flags_edges[cell] = flag;
    } else {
        unreachable();
    }
}





/* Count number of elements */

int MeshSimplicial1D::count_edges() const
{
    return counter_edges;
}

int MeshSimplicial1D::count_vertices() const
{
    return counter_vertices;
}



/* subsimplex relation of edges and vertices */

bool MeshSimplicial1D::contains_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    
    return ( data_edge_vertices[e][0] == v ) || ( data_edge_vertices[e][1] == v );
} 

int MeshSimplicial1D::indexof_edge_vertex( int e, int v ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= v && v < counter_vertices );
    if     ( data_edge_vertices[e][0] == v ) return 0;
    else if( data_edge_vertices[e][1] == v ) return 1;
    else                                     unreachable();
} 

int MeshSimplicial1D::get_edge_vertex( int e, int vi ) const
{
    assert( 0 <= e && e < counter_edges );
    assert( 0 <= vi && vi < 2 );
    return data_edge_vertices[e][vi];
}

const std::array<int,2> MeshSimplicial1D::get_edge_vertices( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    return data_edge_vertices[e];
} 




/* edge parents of a vertex */

int MeshSimplicial1D::count_vertex_edge_parents( int v ) const
{
  return get_edge_parents_of_vertex( v ).size();
}

int MeshSimplicial1D::get_vertex_firstparent_edge( int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  return data_vertex_firstparent_edge[ v ];
}

int MeshSimplicial1D::get_vertex_nextparent_edge( int v, int e ) const
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

int MeshSimplicial1D::get_edge_nextparent_of_vertex( int e, int vi ) const
{
  assert( 0 <= e  && e  < counter_edges    );
  assert( 0 <= vi && vi < 2 );
  return data_edge_nextparents_of_vertices[e][vi];
}


bool MeshSimplicial1D::is_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < counter_vertices );
  assert( 0 <= e && e < counter_edges    );
  return data_edge_vertices[e][0] == v || data_edge_vertices[e][1] == v;
}

int MeshSimplicial1D::indexof_edge_vertex_parent( int e, int v ) const
{
  assert( 0 <= v && v < count_vertices() );
  std::vector<int> edges = get_edge_parents_of_vertex( v );
  
  auto iter = std::find( edges.begin(), edges.end(), e ); 
  assert( iter != edges.end() );
  
  return iter - edges.begin();
}

std::vector<int> MeshSimplicial1D::get_edge_parents_of_vertex( int v ) const
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



void MeshSimplicial1D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    // check();
    
    /* Collect the old data */
    
    int vertex_back       = data_edge_vertices[e][0];
    int vertex_front      = data_edge_vertices[e][1];
    int nextparent_back   = data_edge_nextparents_of_vertices[e][0];
    int nextparent_front  = data_edge_nextparents_of_vertices[e][1];
    int firstparent_back  = data_vertex_firstparent_edge[vertex_back ];
    int firstparent_front = data_vertex_firstparent_edge[vertex_front];
    
    int back_previousparent            = nullindex;
    int back_previousparent_localindex = nullindex;
    if( e != get_vertex_firstparent_edge( vertex_back ) )
      for( back_previousparent = get_vertex_firstparent_edge( vertex_back  );
           back_previousparent != nullindex && get_vertex_nextparent_edge( vertex_back, back_previousparent ) != e; 
           back_previousparent = get_vertex_nextparent_edge( vertex_back, back_previousparent) 
         ); 
    if( back_previousparent  != nullindex ) 
      back_previousparent_localindex  = indexof_edge_vertex( back_previousparent,  vertex_back  );
    
    int front_previousparent  = nullindex;
    int front_previousparent_localindex = nullindex;
    if( e != get_vertex_firstparent_edge( vertex_front ) )
      for( front_previousparent = get_vertex_firstparent_edge( vertex_front );
           front_previousparent != nullindex && get_vertex_nextparent_edge( vertex_front, front_previousparent ) != e; 
           front_previousparent = get_vertex_nextparent_edge( vertex_front, front_previousparent ) 
         ); 
    if( front_previousparent != nullindex )
      front_previousparent_localindex = indexof_edge_vertex( front_previousparent, vertex_front );
    
    /* Assemble the data */
    
    FloatVector midcoordinate = get_edge_midpoint( e );
    
    int ne = counter_edges;
    int nv = counter_vertices;
    
    int back_backvertex       = vertex_back;
    int back_frontvertex      = nv;
    int front_backvertex      = nv;
    int front_frontvertex     = vertex_front;
    
    int back_backnextparent   = nextparent_back;
    int back_frontnextparent  = ne;
    int front_backnextparent  = nullindex;
    int front_frontnextparent = nextparent_front;
    
    int firstparent_newvertex = nv;
    
    /* Allocate memory */
    
    data_edge_nextparents_of_vertices.resize( counter_edges    + 1 );
    data_edge_vertices.resize               ( counter_edges    + 1 );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1 );
    
    flags_edges.resize   ( counter_edges    + 1, SimplexFlag::SimplexFlagInvalid );
    flags_vertices.resize( counter_vertices + 1, SimplexFlag::SimplexFlagInvalid );
    
    /* Write in the data */
    
    data_edge_vertices[e ][0] = back_backvertex;
    data_edge_vertices[e ][1] = back_frontvertex;
    data_edge_vertices[ne][0] = front_backvertex;
    data_edge_vertices[ne][1] = front_frontvertex;
    
    data_edge_nextparents_of_vertices[e ][0] = back_backnextparent;
    data_edge_nextparents_of_vertices[e ][1] = back_frontnextparent;
    data_edge_nextparents_of_vertices[ne][0] = front_backnextparent;
    data_edge_nextparents_of_vertices[ne][1] = front_frontnextparent;
    
    data_vertex_firstparent_edge[nv] = e;
    
    if( back_previousparent  != nullindex ) {
      
      assert( data_edge_nextparents_of_vertices[ back_previousparent ][ back_previousparent_localindex ] == e );
      data_edge_nextparents_of_vertices[ back_previousparent ][ back_previousparent_localindex ] = e;
      
    } else {
      
      assert( data_vertex_firstparent_edge[ vertex_back ] == e );
      data_vertex_firstparent_edge[ vertex_back ] = e;
      
    }
    
    if( front_previousparent != nullindex ) {
      
      assert( data_edge_nextparents_of_vertices[ front_previousparent ][ front_previousparent_localindex ] == ne );
      data_edge_nextparents_of_vertices[ front_previousparent ][ front_previousparent_localindex ] = ne;
      
    } else {
      
      assert( data_vertex_firstparent_edge[ vertex_front ] == e );
      data_vertex_firstparent_edge[ vertex_front ] = ne;
      
    }    
    
    getcoordinates().append( midcoordinate );

    
    flags_edges   [ne] = flags_edges[e];
    flags_vertices[nv] = flags_edges[e];
    
    
    /* Update counter */
    counter_edges++;
    counter_vertices++;
    
    /* Done */
    
    
    // check();
    
}


void MeshSimplicial1D::uniformrefinement()
{
    int old_counter_edges    = counter_edges;
    int old_counter_vertices = counter_vertices;
    
    check();
    
    data_edge_nextparents_of_vertices.reserve( 2 * old_counter_edges );
    data_edge_vertices.reserve               ( 2 * old_counter_edges );
    data_vertex_firstparent_edge.reserve     ( old_counter_vertices  );
    
    getcoordinates().addcapacity( old_counter_edges     );
    
    for( int e = 0; e < old_counter_edges; e++ )
      bisect_edge( e );
      
    check();
}



void MeshSimplicial1D::improved_uniformrefinement()
{
    check();
    
    /* resize the arrays */
    
    data_edge_nextparents_of_vertices.resize( counter_edges * 2                );
    data_edge_vertices.resize               ( counter_edges * 2                );
    data_vertex_firstparent_edge.resize     ( counter_edges + counter_vertices );
    
    getcoordinates().addcoordinates         ( counter_edges                    );
    
    flags_edges.resize   ( counter_edges    * 2,             SimplexFlag::SimplexFlagInvalid );
    flags_vertices.resize( counter_vertices + counter_edges, SimplexFlag::SimplexFlagInvalid );
    
    /* create the new coordinates and fill them up */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      getcoordinates().loadvector( counter_vertices + e, get_edge_midpoint( e ) );
    }
    
    
    /* for each old vertex, set the new parent edge */
    
    for( int v = 0; v < counter_vertices; v++ )
    {
      int p = data_vertex_firstparent_edge[v];
      
      assert( p != nullindex );
      
      int vi = data_edge_vertices[p][0] == v ? 0 : 1;
      
      assert( data_edge_vertices[p][0] == v || data_edge_vertices[p][1] == v );
      assert( data_edge_vertices[p][vi] == v );
      
      data_vertex_firstparent_edge[v] = p + vi * counter_edges;
    }
    
    
    /* for each old edge, relocate the data of the old vertices' old parent edges */
    
    for( int e  = 0; e  < counter_edges;  e++ )
    for( int vi = 0; vi <             2; vi++ )
    {
      int q = data_edge_nextparents_of_vertices[e][vi];
      
      int v = data_edge_vertices[e][vi];
      
      if( q == nullindex ) {
        
        data_edge_nextparents_of_vertices[ e + vi * counter_edges ][vi] = nullindex;
        
      } else {
        
        assert( 0 <= q && q < counter_edges );
        
        int vinp = data_edge_vertices[q][0] == v ? 0 : 1;
        
        assert( data_edge_vertices[q][0] == v || data_edge_vertices[q][1] == v );
        assert( data_edge_vertices[q][vinp] == v );
        
        data_edge_nextparents_of_vertices[ e + vi * counter_edges ][vi] = q + vinp * counter_edges;
      
      } 
      
    }
    
    
    
    /* for each new vertex, set the first and second parent edge from the old edges */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      data_vertex_firstparent_edge[counter_vertices + e] = e;
      
      data_edge_nextparents_of_vertices[e                ][1] = e + counter_edges;
      data_edge_nextparents_of_vertices[e + counter_edges][0] = nullindex;
    }
    
    
    /* for each edge, set the vertices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      int vertex_back  = data_edge_vertices[e][0];
      int vertex_front = data_edge_vertices[e][1];
      
      data_edge_vertices[e                ][0] = vertex_back;
      data_edge_vertices[e                ][1] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][0] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][1] = vertex_front;
    }
    
    
    
    /* set the flags of the newly created simplices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
        flags_edges   [ counter_edges    + e ] = flags_edges[e];
        flags_vertices[ counter_vertices + e ] = flags_edges[e];
    }
    
    
    /* update the counters */
    
    counter_vertices += counter_edges;
    counter_edges += counter_edges;
    
    
    /* DONE */
    
    check();
}



void MeshSimplicial1D::midpoint_refinement( int e )
{
    // check();
    
    bisect_edge( e );
    
    // check();
}

void MeshSimplicial1D::midpoint_refinement_global()
{
    check();
    
    improved_uniformrefinement();
    
    check();
}






FloatVector MeshSimplicial1D::get_edge_midpoint    ( int e ) const
{
    assert( 0 <= e && e < counter_edges );
    FloatVector mid( getouterdimension() );
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] = (   getcoordinates().getdata( get_edge_vertices(e)[0], d ) 
                 + getcoordinates().getdata( get_edge_vertices(e)[1], d )
                ) / 2.;
    return mid;
}





void MeshSimplicial1D::merge( const MeshSimplicial1D& mesh )
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
    
    
    // counters 
    
    counter_vertices += mesh.counter_vertices;
    counter_edges    += mesh.counter_edges;
    
    getcoordinates().append( mesh.getcoordinates() );
    
    check();
}








std::size_t MeshSimplicial1D::memorysize() const
{
    std::size_t ret = 0;

    ret += getcoordinates().memorysize();

    ret += sizeof(getinnerdimension());
    ret += sizeof(getouterdimension());

    ret += sizeof( counter_edges    );
    ret += sizeof( counter_vertices );

    { const auto& D = data_edge_vertices;                ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = data_vertex_firstparent_edge;      ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = data_edge_nextparents_of_vertices; ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    
    { const auto& D = flags_edges;    ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };
    { const auto& D = flags_vertices; ret += sizeof(D) + D.size() * sizeof( std::remove_reference<decltype(D)>::type::value_type); };

    return ret;
}













        



