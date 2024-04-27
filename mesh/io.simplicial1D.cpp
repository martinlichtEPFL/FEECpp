
#include <algorithm>
#include <fstream>
#include <istream>
#include <map>
#include <ostream>
#include <string>
#include <vector>
#include <utility>



#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "mesh.simplicial1D.hpp"
#include "io.simplicial1D.hpp"
#include "coordinates.hpp"
#include "io.coordinates.hpp"




void writeMeshSimplicial1D( const char* filename, const MeshSimplicial1D& mesh, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeMeshSimplicial1D( myfile, mesh, sugar );
    myfile.close();
}

MeshSimplicial1D readMeshSimplicial1D( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    MeshSimplicial1D mesh = readMeshSimplicial1D( myfile );
    myfile.close();
    return mesh;
}



void writeMeshSimplicial1D( std::ostream& out, const MeshSimplicial1D& mesh, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Writing simplicial 1D Mesh..." << nl;
    
    if( sugar ) out << "number of edges: " << nl;;
    out << mesh.count_edges() << nl;
    
    if( sugar ) out << "number of vertices: " << nl;
    out << mesh.count_vertices() << nl;
    
    if( sugar ) out << "external dimension: " << nl;
    out << mesh.getcoordinates().getdimension() << nl;
    
    /* edge -> vertices */
    if( sugar ) out << "for each edge, the vertices: " << nl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_vertex( e, 0 )
            << space
            << mesh.get_edge_vertex( e, 1 ) << nl;
    }
    
    if( sugar ) out << "for each vertex, the first parent edge: " << nl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_edge( v )
            << nl;
    }
    
    /* edge -> next parent of vertex vertices */
    if( sugar ) out << "for each edge, the next neighbors: " << nl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_nextparent_of_vertex( e, 0 )
            << space
            << mesh.get_edge_nextparent_of_vertex( e, 1 )
            << nl;
    }
    
    assert( out.good() );

    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << static_cast<int>( mesh.get_flag(1,e) ) << nl;
    }

    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << static_cast<int>( mesh.get_flag(0,v) ) << nl;
    }
    
    assert( out.good() );
    
    writeCoordinates( out, mesh.getcoordinates(), sugar );
}




MeshSimplicial1D readMeshSimplicial1D( std::istream& in )
{
    int counter_edges, counter_vertices, dim;
    
    in >> counter_edges
       >> counter_vertices
       >> dim;
    
    int nullindex = MeshSimplicial1D::nullindex;
    
    /* edge -> vertices */
    
    std::vector< std::array<int,2> > edge_vertices( counter_edges, { nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_edge( counter_vertices, nullindex );
    std::vector< std::array<int,2> > edge_nextparents_of_vertices( counter_edges, { nullindex, nullindex } );
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_vertices[e][0] >> edge_vertices[e][1];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_edge[v];
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_nextparents_of_vertices[e][0] >> edge_nextparents_of_vertices[e][1];
    
    /* flags */

    std::vector<SimplexFlag> flags_edges   ( counter_edges,    SimplexFlag::SimplexFlagInvalid );
    std::vector<SimplexFlag> flags_vertices( counter_vertices, SimplexFlag::SimplexFlagInvalid );
    
    for( int e = 0; e < counter_edges; e++ ) {
        int temp; 
        in >> temp;
        flags_edges[e] = static_cast<SimplexFlag>( temp );
    }

    for( int v = 0; v < counter_vertices; v++ ) {
        int temp; 
        in >> temp;
        flags_vertices[v] = static_cast<SimplexFlag>( temp );
    }
    
    assert( in.good() );
    
    /* coordinates */
    
    Coordinates coords = readCoordinates( in );
    
    assert( in.good() );
    
    /* return */
    
    auto ret = MeshSimplicial1D( dim, coords, edge_vertices, edge_nextparents_of_vertices, vertex_firstparent_edge );

    ret.set_flags( 1, flags_edges    );
    ret.set_flags( 0, flags_vertices );

    return ret;
}



        
        
        
        
        
