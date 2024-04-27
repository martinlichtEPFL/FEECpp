
#include <algorithm>
#include <fstream>
#include <istream>
#include <ostream>
#include <map>
#include <string>
#include <utility>
#include <vector>



#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "mesh.simplicial3D.hpp"
#include "io.simplicial3D.hpp"
#include "coordinates.hpp"
#include "io.coordinates.hpp"




void writeMeshSimplicial3D( const char* filename, const MeshSimplicial3D& mesh, bool sugar )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    writeMeshSimplicial3D( myfile, mesh, sugar );
    myfile.close();
}

MeshSimplicial3D readMeshSimplicial3D( const char* filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    MeshSimplicial3D mesh = readMeshSimplicial3D( myfile );
    myfile.close();
    return mesh;
}



void writeMeshSimplicial3D( std::ostream& out, const MeshSimplicial3D& mesh, bool sugar )
{
    /* Preamble */
    if( sugar ) out << "Writing simplicial 3D Mesh..." << nl;
    
    if( sugar ) out << "number of tetrahedra: " << nl;;
    out << mesh.count_tetrahedra() << nl;
    
    if( sugar ) out << "number of faces: " << nl;;
    out << mesh.count_faces() << nl;
    
    if( sugar ) out << "number of edges: " << nl;;
    out << mesh.count_edges() << nl;
    
    if( sugar ) out << "number of vertices: " << nl;
    out << mesh.count_vertices() << nl;
    
    if( sugar ) out << "external dimension: " << nl;
    out << mesh.getcoordinates().getdimension() << nl;
    
    /* tetrahedron -> faces */
    
    if( sugar ) out << "for each tetrahedron, the faces: " << nl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_face( t, 0 )
            << space
            << mesh.get_tetrahedron_face( t, 1 )
            << space
            << mesh.get_tetrahedron_face( t, 2 )
            << space
            << mesh.get_tetrahedron_face( t, 3 ) << nl;
    }
    
    if( sugar ) out << "for each face, the first parent tetrahedron: " << nl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_firstparent_tetrahedron( f )
            << nl;
    }
    
    if( sugar ) out << "for each tetrahedron, the next neighbors: " << nl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_nextparent_of_face( t, 0 )
            << space
            << mesh.get_tetrahedron_nextparent_of_face( t, 1 )
            << space
            << mesh.get_tetrahedron_nextparent_of_face( t, 2 )
            << space
            << mesh.get_tetrahedron_nextparent_of_face( t, 3 )
            << nl;
    }
    
    /* tetrahedron -> edges */
    
    if( sugar ) out << "for each tetrahedron, the edges: " << nl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_edge( t, 0 )
            << space
            << mesh.get_tetrahedron_edge( t, 1 )
            << space
            << mesh.get_tetrahedron_edge( t, 2 )
            << space
            << mesh.get_tetrahedron_edge( t, 3 )
            << space
            << mesh.get_tetrahedron_edge( t, 4 )
            << space
            << mesh.get_tetrahedron_edge( t, 5 ) << nl;
    }
    
    if( sugar ) out << "for each edge, the first parent tetrahedron: " << nl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_firstparent_tetrahedron( e )
            << nl;
    }
    
    if( sugar ) out << "for each tetrahedron, the next neighbors: " << nl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_nextparent_of_edge( t, 0 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 1 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 2 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 3 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 4 )
            << space
            << mesh.get_tetrahedron_nextparent_of_edge( t, 5 )
            << nl;
    }
    
    /* tetrahedron -> vertices */
    
    if( sugar ) out << "for each tetrahedron, the vertices: " << nl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_vertex( t, 0 )
            << space
            << mesh.get_tetrahedron_vertex( t, 1 )
            << space
            << mesh.get_tetrahedron_vertex( t, 2 )
            << space
            << mesh.get_tetrahedron_vertex( t, 3 ) << nl;
    }
    
    if( sugar ) out << "for each vertex, the first parent tetrahedron: " << nl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_tetrahedron( v )
            << nl;
    }
    
    if( sugar ) out << "for each tetrahedron, the next neighbors: " << nl;
    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << mesh.get_tetrahedron_nextparent_of_vertex( t, 0 )
            << space
            << mesh.get_tetrahedron_nextparent_of_vertex( t, 1 )
            << space
            << mesh.get_tetrahedron_nextparent_of_vertex( t, 2 )
            << space
            << mesh.get_tetrahedron_nextparent_of_vertex( t, 3 )
            << nl;
    }
    
    assert( out.good() );
    
    /* face -> edges */
    
    if( sugar ) out << "for each face, the edges: " << nl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_edge( f, 0 )
            << space
            << mesh.get_face_edge( f, 1 )
            << space
            << mesh.get_face_edge( f, 2 ) << nl;
    }
    
    if( sugar ) out << "for each edge, the first parent face: " << nl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_firstparent_face( e )
            << nl;
    }
    
    if( sugar ) out << "for each face, the next neighbors: " << nl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_nextparent_of_edge( f, 0 )
            << space
            << mesh.get_face_nextparent_of_edge( f, 1 )
            << space
            << mesh.get_face_nextparent_of_edge( f, 2 )
            << nl;
    }
    
    /* face -> vertices */
    
    if( sugar ) out << "for each face, the vertices: " << nl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_vertex( f, 0 )
            << space
            << mesh.get_face_vertex( f, 1 )
            << space
            << mesh.get_face_vertex( f, 2 ) << nl;
    }
    
    if( sugar ) out << "for each vertex, the first parent face: " << nl;
    for( int v = 0; v < mesh.count_vertices(); v++ ) {
        if( sugar ) out << v << ": ";
        out << mesh.get_vertex_firstparent_face( v )
            << nl;
    }
    
    if( sugar ) out << "for each face, the next neighbors: " << nl;
    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << mesh.get_face_nextparent_of_vertex( f, 0 )
            << space
            << mesh.get_face_nextparent_of_vertex( f, 1 )
            << space
            << mesh.get_face_nextparent_of_vertex( f, 2 )
            << nl;
    }
    
    assert( out.good() );
    
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
    
    if( sugar ) out << "for each edge, the next neighbors: " << nl;
    for( int e = 0; e < mesh.count_edges(); e++ ) {
        if( sugar ) out << e << ": ";
        out << mesh.get_edge_nextparent_of_vertex( e, 0 )
            << space
            << mesh.get_edge_nextparent_of_vertex( e, 1 )
            << nl;
    }
    
    assert( out.good() );

    for( int t = 0; t < mesh.count_tetrahedra(); t++ ) {
        if( sugar ) out << t << ": ";
        out << static_cast<int>( mesh.get_flag(3,t) ) << nl;
    }

    for( int f = 0; f < mesh.count_faces(); f++ ) {
        if( sugar ) out << f << ": ";
        out << static_cast<int>( mesh.get_flag(2,f) ) << nl;
    }

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





MeshSimplicial3D readMeshSimplicial3D( std::istream& in )
{
    int counter_tetrahedra, counter_faces, counter_edges, counter_vertices, dim;
    
    in >> counter_tetrahedra
       >> counter_faces
       >> counter_edges
       >> counter_vertices
       >> dim;
    
    int nullindex = MeshSimplicial3D::nullindex;
    
    /* tetrahedron -> faces */
    
    std::vector< std::array<int,4> > tetrahedron_faces( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    std::vector< int               > face_firstparent_tetrahedron( counter_faces, nullindex );
    std::vector< std::array<int,4> > tetrahedron_nextparents_of_faces( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_faces[t][0]
           >> tetrahedron_faces[t][1]
           >> tetrahedron_faces[t][2]
           >> tetrahedron_faces[t][3];
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_firstparent_tetrahedron[f];
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_nextparents_of_faces[t][0] 
           >> tetrahedron_nextparents_of_faces[t][1] 
           >> tetrahedron_nextparents_of_faces[t][2] 
           >> tetrahedron_nextparents_of_faces[t][3];
    
    assert( in.good() );
    
    /* tetrahedron -> edges */
    
    std::vector< std::array<int,6> > tetrahedron_edges( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    std::vector< int               > edge_firstparent_tetrahedron( counter_edges, nullindex );
    std::vector< std::array<int,6> > tetrahedron_nextparents_of_edges( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_edges[t][0]
           >> tetrahedron_edges[t][1]
           >> tetrahedron_edges[t][2]
           >> tetrahedron_edges[t][3]
           >> tetrahedron_edges[t][4]
           >> tetrahedron_edges[t][5];
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_firstparent_tetrahedron[e];
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_nextparents_of_edges[t][0] 
           >> tetrahedron_nextparents_of_edges[t][1] 
           >> tetrahedron_nextparents_of_edges[t][2] 
           >> tetrahedron_nextparents_of_edges[t][3] 
           >> tetrahedron_nextparents_of_edges[t][4] 
           >> tetrahedron_nextparents_of_edges[t][5];
    
    assert( in.good() );
    
    /* tetrahedron -> vertices */
    
    std::vector< std::array<int,4> > tetrahedron_vertices( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_tetrahedron( counter_vertices, nullindex );
    std::vector< std::array<int,4> > tetrahedron_nextparents_of_vertices( counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_vertices[t][0]
           >> tetrahedron_vertices[t][1]
           >> tetrahedron_vertices[t][2]
           >> tetrahedron_vertices[t][3];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_tetrahedron[v];
    
    for( int t = 0; t < counter_tetrahedra; t++ )
        in >> tetrahedron_nextparents_of_vertices[t][0] 
           >> tetrahedron_nextparents_of_vertices[t][1] 
           >> tetrahedron_nextparents_of_vertices[t][2] 
           >> tetrahedron_nextparents_of_vertices[t][3];
    
    assert( in.good() );
    
    /* face -> edges */
    
    std::vector< std::array<int,3> > face_edges( counter_faces, { nullindex, nullindex, nullindex } );
    std::vector< int               > edge_firstparent_face( counter_edges, nullindex );
    std::vector< std::array<int,3> > face_nextparents_of_edges( counter_faces, { nullindex, nullindex, nullindex } );
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_edges[f][0]
           >> face_edges[f][1]
           >> face_edges[f][2];
    
    for( int v = 0; v < counter_edges; v++ )
        in >> edge_firstparent_face[v];
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_nextparents_of_edges[f][0] 
           >> face_nextparents_of_edges[f][1] 
           >> face_nextparents_of_edges[f][2];
    
    assert( in.good() );
    
    /* face -> vertices */
    
    std::vector< std::array<int,3> > face_vertices( counter_faces, { nullindex, nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_face( counter_vertices, nullindex );
    std::vector< std::array<int,3> > face_nextparents_of_vertices( counter_faces, { nullindex, nullindex, nullindex } );
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_vertices[f][0]
           >> face_vertices[f][1]
           >> face_vertices[f][2];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_face[v];
    
    for( int f = 0; f < counter_faces; f++ )
        in >> face_nextparents_of_vertices[f][0] 
           >> face_nextparents_of_vertices[f][1] 
           >> face_nextparents_of_vertices[f][2];
    
    assert( in.good() );
    
    /* edges -> vertices */
    
    std::vector< std::array<int,2> > edge_vertices( counter_edges, { nullindex, nullindex } );
    std::vector< int               > vertex_firstparent_edge( counter_vertices, nullindex );
    std::vector< std::array<int,2> > edge_nextparents_of_vertices( counter_edges, { nullindex, nullindex } );
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_vertices[e][0] >> edge_vertices[e][1];
    
    for( int v = 0; v < counter_vertices; v++ )
        in >> vertex_firstparent_edge[v];
    
    for( int e = 0; e < counter_edges; e++ )
        in >> edge_nextparents_of_vertices[e][0] >> edge_nextparents_of_vertices[e][1];
    
    assert( in.good() );
    
    /* flags */

    std::vector<SimplexFlag> flags_tetrahedra( counter_tetrahedra, SimplexFlag::SimplexFlagInvalid );
    std::vector<SimplexFlag> flags_faces     ( counter_faces,      SimplexFlag::SimplexFlagInvalid );
    std::vector<SimplexFlag> flags_edges     ( counter_edges,      SimplexFlag::SimplexFlagInvalid );
    std::vector<SimplexFlag> flags_vertices  ( counter_vertices,   SimplexFlag::SimplexFlagInvalid );
    
    for( int t = 0; t < counter_tetrahedra; t++ ) {
        int temp; 
        in >> temp;
        flags_tetrahedra[t] = static_cast<SimplexFlag>( temp );
    }

    for( int f = 0; f < counter_faces; f++ ) {
        int temp; 
        in >> temp;
        flags_faces[f] = static_cast<SimplexFlag>( temp );
    }

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
    
    auto ret = MeshSimplicial3D( dim, coords, 
                             tetrahedron_faces,    face_firstparent_tetrahedron,   tetrahedron_nextparents_of_faces, 
                             tetrahedron_edges,    edge_firstparent_tetrahedron,   tetrahedron_nextparents_of_edges, 
                             tetrahedron_vertices, vertex_firstparent_tetrahedron, tetrahedron_nextparents_of_vertices, 
                             face_edges,           edge_firstparent_face,          face_nextparents_of_edges, 
                             face_vertices,        vertex_firstparent_face,        face_nextparents_of_vertices, 
                             edge_vertices,        vertex_firstparent_edge,        edge_nextparents_of_vertices
                           );

    ret.set_flags( 3, flags_tetrahedra );
    ret.set_flags( 2, flags_faces      );
    ret.set_flags( 1, flags_edges      );
    ret.set_flags( 0, flags_vertices   );

    return ret;
}



