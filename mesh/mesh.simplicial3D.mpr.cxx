
#include <string>
#include <vector>
#include <map>
#include <utility>
#include <algorithm>
#include <iostream>
#include <fstream>


#include "../basic.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "coordinates.hpp"
#include "mesh.hpp"

#include "mesh.simplicial3D.hpp"










void MeshSimplicial3D::midpoint_refinement( int t )
{
    // check();
    
//     assert( 0 <= t && t < counter_faces );
//     
//     /* Allocate memory */
//     
//     data_tetrahedron_faces.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
//     data_face_firstparent_tetrahedron.resize    ( counter_faces + 6,                                         nullindex   );
//     data_tetrahedron_nextparents_of_faces.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
//     
//     data_tetrahedron_edges.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
//     data_edge_firstparent_tetrahedron.resize    ( counter_edges + 4,                                                               nullindex   );
//     data_tetrahedron_nextparents_of_edges.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
//     
//     data_tetrahedron_vertices.resize               ( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
//     data_vertex_firstparent_tetrahedron.resize     ( counter_vertices + 1,                                      nullindex   );
//     data_tetrahedron_nextparents_of_vertices.resize( counter_tetrahedra + 3, { nullindex, nullindex, nullindex, nullindex } );
//     
//     data_face_edges.resize               ( counter_faces + 6, { nullindex, nullindex, nullindex } );
//     data_edge_firstparent_face.resize    ( counter_edges + 4,                         nullindex   );
//     data_face_nextparents_of_edges.resize( counter_faces + 6, { nullindex, nullindex, nullindex } );
//     
//     data_face_vertices.resize               ( counter_faces + 6, { nullindex, nullindex, nullindex } );
//     data_vertex_firstparent_face.resize     ( counter_vertices + 1,                      nullindex   );
//     data_face_nextparents_of_vertices.resize( counter_faces + 6, { nullindex, nullindex, nullindex } );
//     
//     data_edge_vertices.resize               ( counter_edges + 4,   { nullindex, nullindex } );
//     data_vertex_firstparent_edge.resize     ( counter_vertices + 1,             nullindex   );
//     data_edge_nextparents_of_vertices.resize( counter_edges + 4,   { nullindex, nullindex } );
//     
//     getcoordinates().addcoordinates( 1 );
//     
//     
//     /* load the new coordinate */
//     
//     getcoordinates().loadvector( counter_vertices, get_face_midpoint( t ) );
//     
//     
//     /* assemble the data and auxiliary variables */
//     int t_f0 = data_tetrahedron_faces[t][0];
//     int t_f1 = data_tetrahedron_faces[t][1];
//     int t_f2 = data_tetrahedron_faces[t][2];
//     int t_f3 = data_tetrahedron_faces[t][3];
//     
//     int t_e0 = data_tetrahedron_edges[t][0];
//     int t_e1 = data_tetrahedron_edges[t][1];
//     int t_e2 = data_tetrahedron_edges[t][2];
//     int t_e3 = data_tetrahedron_edges[t][3];
//     int t_e4 = data_tetrahedron_edges[t][4];
//     int t_e5 = data_tetrahedron_edges[t][5];
//     
//     int t_v0 = data_tetrahedron_vertices[t][0];
//     int t_v1 = data_tetrahedron_vertices[t][1];
//     int t_v2 = data_tetrahedron_vertices[t][2];
//     int t_v3 = data_tetrahedron_vertices[t][3];
//     
//     int t_f_n0 = data_tetrahedron_nextparents_of_faces[t][0];
//     int t_f_n1 = data_tetrahedron_nextparents_of_faces[t][1];
//     int t_f_n2 = data_tetrahedron_nextparents_of_faces[t][2];
//     int t_f_n3 = data_tetrahedron_nextparents_of_faces[t][3];
//     
//     int t_e_n0 = data_tetrahedron_nextparents_of_edges[t][0];
//     int t_e_n1 = data_tetrahedron_nextparents_of_edges[t][1];
//     int t_e_n2 = data_tetrahedron_nextparents_of_edges[t][2];
//     int t_e_n3 = data_tetrahedron_nextparents_of_edges[t][3];
//     int t_e_n4 = data_tetrahedron_nextparents_of_edges[t][4];
//     int t_e_n5 = data_tetrahedron_nextparents_of_edges[t][5];
//     
//     int t_v_n0 = data_tetrahedron_nextparents_of_vertices[t][0];
//     int t_v_n1 = data_tetrahedron_nextparents_of_vertices[t][1];
//     int t_v_n2 = data_tetrahedron_nextparents_of_vertices[t][2];
//     int t_v_n3 = data_tetrahedron_nextparents_of_vertices[t][3];
//     
//     
//     int t0 = t;
//     int t1 = counter_faces;
//     int t2 = counter_faces + 1;
//     int t3 = counter_faces + 2;
//     
//     int f0 = counter_faces + 0;
//     int f1 = counter_faces + 1;
//     int f2 = counter_faces + 2;
//     int f3 = counter_faces + 3;
//     int f4 = counter_faces + 4;
//     int f5 = counter_faces + 5;
//     
//     int e0 = counter_edges + 0;
//     int e1 = counter_edges + 1;
//     int e2 = counter_edges + 2;
//     int e3 = counter_edges + 3;
//     
//     int vn = counter_vertices;
//     
//     
//     /* fill in: downward */ // TODO
//     data_tetrahedron_vertices[ t0 ][0] = vn;
//     data_tetrahedron_vertices[ t0 ][1] = t_v0;
//     data_tetrahedron_vertices[ t0 ][2] = t_v1;
//     data_tetrahedron_vertices[ t0 ][3] = t_v2;
//     
//     data_tetrahedron_vertices[ t1 ][0] = vn;
//     data_tetrahedron_vertices[ t1 ][1] = t_v0;
//     data_tetrahedron_vertices[ t1 ][2] = t_v1;
//     data_tetrahedron_vertices[ t1 ][3] = t_v3;
//     
//     data_tetrahedron_vertices[ t2 ][0] = vn;
//     data_tetrahedron_vertices[ t2 ][1] = t_v0;
//     data_tetrahedron_vertices[ t2 ][2] = t_v2;
//     data_tetrahedron_vertices[ t2 ][3] = t_v3;
//     
//     data_tetrahedron_vertices[ t3 ][0] = vn;
//     data_tetrahedron_vertices[ t3 ][1] = t_v1;
//     data_tetrahedron_vertices[ t3 ][2] = t_v2;
//     data_tetrahedron_vertices[ t3 ][3] = t_v3;
//     
//     data_tetrahedron_edges[ t0 ][0] = e0;
//     data_tetrahedron_edges[ t0 ][1] = e1;
//     data_tetrahedron_edges[ t0 ][2] = e2;
//     data_tetrahedron_edges[ t0 ][3] = t_e0;
//     data_tetrahedron_edges[ t0 ][4] = t_e1;
//     data_tetrahedron_edges[ t0 ][5] = t_e3;
//     
//     data_tetrahedron_edges[ t1 ][0] = e0;
//     data_tetrahedron_edges[ t1 ][1] = e1;
//     data_tetrahedron_edges[ t1 ][2] = e3;
//     data_tetrahedron_edges[ t1 ][3] = t_e0;
//     data_tetrahedron_edges[ t1 ][4] = t_e2;
//     data_tetrahedron_edges[ t1 ][5] = t_e4;
//     
//     data_tetrahedron_edges[ t2 ][0] = e0;
//     data_tetrahedron_edges[ t2 ][1] = e2;
//     data_tetrahedron_edges[ t2 ][2] = e3;
//     data_tetrahedron_edges[ t2 ][3] = t_e1;
//     data_tetrahedron_edges[ t2 ][4] = t_e2;
//     data_tetrahedron_edges[ t2 ][5] = t_e5;
//     
//     data_tetrahedron_edges[ t3 ][0] = e1;
//     data_tetrahedron_edges[ t3 ][1] = e2;
//     data_tetrahedron_edges[ t3 ][2] = e3;
//     data_tetrahedron_edges[ t3 ][3] = t_e3;
//     data_tetrahedron_edges[ t3 ][4] = t_e4;
//     data_tetrahedron_edges[ t3 ][5] = t_e5;
//     
//     data_tetrahedron_faces[ t0 ][0] = f0;
//     data_tetrahedron_faces[ t0 ][1] = f1;
//     data_tetrahedron_faces[ t0 ][2] = f3;
//     data_tetrahedron_faces[ t0 ][3] = t_f0;
//     
//     data_tetrahedron_faces[ t1 ][0] = f0;
//     data_tetrahedron_faces[ t1 ][1] = f2;
//     data_tetrahedron_faces[ t1 ][2] = f4;
//     data_tetrahedron_faces[ t1 ][3] = t_f1;
//     
//     data_tetrahedron_faces[ t2 ][0] = f1;
//     data_tetrahedron_faces[ t2 ][1] = f2;
//     data_tetrahedron_faces[ t2 ][2] = f5;
//     data_tetrahedron_faces[ t2 ][3] = t_f2;
//     
//     data_tetrahedron_faces[ t3 ][0] = f3;
//     data_tetrahedron_faces[ t3 ][1] = f4;
//     data_tetrahedron_faces[ t3 ][2] = f5;
//     data_tetrahedron_faces[ t3 ][3] = t_f3;
//     
//     data_face_vertices[ f0 ][0] = vn;
//     data_face_vertices[ f0 ][1] = t_v0;
//     data_face_vertices[ f0 ][2] = t_v1;
//     
//     data_face_vertices[ f1 ][0] = vn;
//     data_face_vertices[ f1 ][1] = t_v0;
//     data_face_vertices[ f1 ][2] = t_v2;
//     
//     data_face_vertices[ f2 ][0] = vn;
//     data_face_vertices[ f2 ][1] = t_v0;
//     data_face_vertices[ f2 ][2] = t_v3;
//     
//     data_face_vertices[ f3 ][0] = vn;
//     data_face_vertices[ f3 ][1] = t_v1;
//     data_face_vertices[ f3 ][2] = t_v2;
//     
//     data_face_vertices[ f4 ][0] = vn;
//     data_face_vertices[ f4 ][1] = t_v1;
//     data_face_vertices[ f4 ][2] = t_v3;
//     
//     data_face_vertices[ f5 ][0] = vn;
//     data_face_vertices[ f5 ][1] = t_v2;
//     data_face_vertices[ f5 ][2] = t_v3;
//     
//     data_face_edges[ f0 ][0] = e0;
//     data_face_edges[ f0 ][1] = e1;
//     data_face_edges[ f0 ][2] = t_e0;
//     
//     data_face_edges[ f1 ][0] = e0;
//     data_face_edges[ f1 ][1] = e2;
//     data_face_edges[ f1 ][2] = t_e1;
//     
//     data_face_edges[ f2 ][0] = e0;
//     data_face_edges[ f2 ][1] = e3;
//     data_face_edges[ f2 ][2] = t_e2;
//     
//     data_face_edges[ f3 ][0] = e1;
//     data_face_edges[ f3 ][1] = e2;
//     data_face_edges[ f3 ][2] = t_e3;
//     
//     data_face_edges[ f4 ][0] = e1;
//     data_face_edges[ f4 ][1] = e3;
//     data_face_edges[ f4 ][2] = t_e4;
//     
//     data_face_edges[ f5 ][0] = e2;
//     data_face_edges[ f5 ][1] = e3;
//     data_face_edges[ f5 ][2] = t_e5;
//     
//     data_edge_vertices[ e0 ][0] = vn;
//     data_edge_vertices[ e0 ][1] = t_v0;
//     
//     data_edge_vertices[ e1 ][0] = vn;
//     data_edge_vertices[ e1 ][1] = t_v1;
//     
//     data_edge_vertices[ e2 ][0] = vn;
//     data_edge_vertices[ e2 ][1] = t_v2;
//     
//     data_edge_vertices[ e3 ][0] = vn;
//     data_edge_vertices[ e3 ][1] = t_v3;
//     
//     /* fill in: parentlist */
//     data_tetrahedron_nextparents_of_vertices[ t0 ][0] = t1;
//     data_tetrahedron_nextparents_of_vertices[ t0 ][1] = t1;
//     data_tetrahedron_nextparents_of_vertices[ t0 ][2] = t1;
//     data_tetrahedron_nextparents_of_vertices[ t0 ][3] = t2;
//     
//     data_tetrahedron_nextparents_of_vertices[ t1 ][0] = t2;
//     data_tetrahedron_nextparents_of_vertices[ t1 ][1] = t2;
//     data_tetrahedron_nextparents_of_vertices[ t1 ][2] = t3;
//     data_tetrahedron_nextparents_of_vertices[ t1 ][3] = t2;
//     
//     data_tetrahedron_nextparents_of_vertices[ t2 ][0] = t3;
//     data_tetrahedron_nextparents_of_vertices[ t2 ][1] = data_vertex_firstparent_tetrahedron[ t_v0 ];
//     data_tetrahedron_nextparents_of_vertices[ t2 ][2] = t3;
//     data_tetrahedron_nextparents_of_vertices[ t2 ][3] = t3;
//     data_vertex_firstparent_tetrahedron[ t_v0 ] = t0;
//     
//     data_tetrahedron_nextparents_of_vertices[ t3 ][0] = nullindex;
//     data_tetrahedron_nextparents_of_vertices[ t3 ][1] = data_vertex_firstparent_tetrahedron[ t_v1 ];
//     data_tetrahedron_nextparents_of_vertices[ t3 ][2] = data_vertex_firstparent_tetrahedron[ t_v2 ];
//     data_tetrahedron_nextparents_of_vertices[ t3 ][3] = data_vertex_firstparent_tetrahedron[ t_v3 ];
//     data_vertex_firstparent_tetrahedron[ t_v1 ] = t0;
//     data_vertex_firstparent_tetrahedron[ t_v2 ] = t0;
//     data_vertex_firstparent_tetrahedron[ t_v3 ] = t1;
//     
//     data_tetrahedron_nextparents_of_edges[ t0 ][0] = t1; // e0
//     data_tetrahedron_nextparents_of_edges[ t0 ][1] = t1; // e1
//     data_tetrahedron_nextparents_of_edges[ t0 ][2] = t2; // e2
//     data_tetrahedron_nextparents_of_edges[ t0 ][3] = t1; // t_e_n0
//     data_tetrahedron_nextparents_of_edges[ t0 ][4] = t2; // t_e_n1
//     data_tetrahedron_nextparents_of_edges[ t0 ][5] = t3 // t_e_n3
//     
//     data_tetrahedron_nextparents_of_edges[ t1 ][0] = t2; // e0
//     data_tetrahedron_nextparents_of_edges[ t1 ][1] = t3; // e1
//     data_tetrahedron_nextparents_of_edges[ t1 ][2] = t2; // e3
//     data_tetrahedron_nextparents_of_edges[ t1 ][3] = t_e_n0; //FIXME: relink below 
//     data_tetrahedron_nextparents_of_edges[ t1 ][4] = t2; // t_e_n2
//     data_tetrahedron_nextparents_of_edges[ t1 ][5] = t3; // t_e_n4
//     
//     data_tetrahedron_nextparents_of_edges[ t2 ][0] = nullindex; data_edge_firstparent_tetrahedron[ e0 ] = t0;
//     data_tetrahedron_nextparents_of_edges[ t2 ][1] = t3; // e2
//     data_tetrahedron_nextparents_of_edges[ t2 ][2] = t3; // e3
//     data_tetrahedron_nextparents_of_edges[ t2 ][3] = t_e_n1; //FIXME: relink below
//     data_tetrahedron_nextparents_of_edges[ t2 ][4] = t_e_n2; //FIXME: relink below
//     data_tetrahedron_nextparents_of_edges[ t2 ][5] = t3; // t_e_n5
//     
//     data_tetrahedron_nextparents_of_edges[ t3 ][0] = nullindex; data_edge_firstparent_tetrahedron[ e1 ] = t0;
//     data_tetrahedron_nextparents_of_edges[ t3 ][1] = nullindex; data_edge_firstparent_tetrahedron[ e2 ] = t0;
//     data_tetrahedron_nextparents_of_edges[ t3 ][2] = nullindex; data_edge_firstparent_tetrahedron[ e3 ] = t0;
//     data_tetrahedron_nextparents_of_edges[ t3 ][3] = t_e_n3; //FIXME: relink below
//     data_tetrahedron_nextparents_of_edges[ t3 ][4] = t_e_n4; //FIXME: relink below
//     data_tetrahedron_nextparents_of_edges[ t3 ][5] = t_e_n5; //FIXME: relink below
//     
//     data_tetrahedron_nextparents_of_faces[ t0 ][0] = t1;
//     data_tetrahedron_nextparents_of_faces[ t0 ][1] = t2;
//     data_tetrahedron_nextparents_of_faces[ t0 ][2] = t3;
//     data_tetrahedron_nextparents_of_faces[ t0 ][3] = t_f_n0;
//     
//     data_tetrahedron_nextparents_of_faces[ t1 ][0] = nullindex;
//     data_tetrahedron_nextparents_of_faces[ t1 ][1] = t2;
//     data_tetrahedron_nextparents_of_faces[ t1 ][2] = t3;
//     data_tetrahedron_nextparents_of_faces[ t1 ][3] = t_f_n1;
//     
//     data_tetrahedron_nextparents_of_faces[ t2 ][0] = nullindex;
//     data_tetrahedron_nextparents_of_faces[ t2 ][1] = nullindex;
//     data_tetrahedron_nextparents_of_faces[ t2 ][2] = t3;
//     data_tetrahedron_nextparents_of_faces[ t2 ][3] = t_f_n2;
//     
//     data_tetrahedron_nextparents_of_faces[ t3 ][0] = nullindex;
//     data_tetrahedron_nextparents_of_faces[ t3 ][1] = nullindex;
//     data_tetrahedron_nextparents_of_faces[ t3 ][2] = nullindex;
//     data_tetrahedron_nextparents_of_faces[ t3 ][3] = t_f_n3;
//     
//     data_face_nextparents_of_vertices[ f0 ][0] = f1;
//     data_face_nextparents_of_vertices[ f0 ][1] = f1; // t_v0
//     data_face_nextparents_of_vertices[ f0 ][2] = f3; // t_v1
//     
//     data_face_nextparents_of_vertices[ f1 ][0] = f2;
//     data_face_nextparents_of_vertices[ f1 ][1] = f2; // t_v0
//     data_face_nextparents_of_vertices[ f1 ][2] = f3; // t_v2
//     
//     data_face_nextparents_of_vertices[ f2 ][0] = f3;
//     data_face_nextparents_of_vertices[ f2 ][1] = data_vertex_firstparent_face[ t_v0 ]; // t_v0
//     data_face_nextparents_of_vertices[ f2 ][2] = f4; // t_v3
//     
//     data_face_nextparents_of_vertices[ f3 ][0] = f4;
//     data_face_nextparents_of_vertices[ f3 ][1] = f4; // t_v1
//     data_face_nextparents_of_vertices[ f3 ][2] = f5; // t_v2
//     
//     data_face_nextparents_of_vertices[ f4 ][0] = f5;
//     data_face_nextparents_of_vertices[ f4 ][1] = data_vertex_firstparent_face[ t_v1 ]; // t_v1
//     data_face_nextparents_of_vertices[ f4 ][2] = f5; // t_v3
//     
//     data_face_nextparents_of_vertices[ f5 ][0] = nullindex; data_vertex_firstparent_face[ vn ] = f0;
//     data_face_nextparents_of_vertices[ f5 ][1] = data_vertex_firstparent_face[ t_v2 ]; // t_v2
//     data_face_nextparents_of_vertices[ f5 ][2] = data_vertex_firstparent_face[ t_v3 ]; // t_v3
//     
//     data_vertex_firstparent_face[ t_v0 ] = f0;
//     data_vertex_firstparent_face[ t_v1 ] = f0;
//     data_vertex_firstparent_face[ t_v2 ] = f1;
//     data_vertex_firstparent_face[ t_v3 ] = f2;
//     
//     data_face_nextparents_of_edges[ f0 ][0] = f1; // e0
//     data_face_nextparents_of_edges[ f0 ][1] = f3; // e1
//     data_face_nextparents_of_edges[ f0 ][2] = data_edge_firstparent_face[ t_e0 ]; // t_e0
//     data_edge_firstparent_face[ t_e0 ] = f0;
//     
//     data_face_nextparents_of_edges[ f1 ][0] = f2; // e0
//     data_face_nextparents_of_edges[ f1 ][1] = f3; // e2
//     data_face_nextparents_of_edges[ f1 ][2] = data_edge_firstparent_face[ t_e1 ]; // t_e1
//     data_edge_firstparent_face[ t_e1 ] = f1;
//     
//     data_face_nextparents_of_edges[ f2 ][0] = nullindex; // e0
//     data_face_nextparents_of_edges[ f2 ][1] = f4; // e3
//     data_face_nextparents_of_edges[ f2 ][2] = data_edge_firstparent_face[ t_e2 ]; // t_e2
//     data_edge_firstparent_face[ t_e2 ] = f2;
//     
//     data_face_nextparents_of_edges[ f3 ][0] = f4; // e1
//     data_face_nextparents_of_edges[ f3 ][1] = f5; // e2
//     data_face_nextparents_of_edges[ f3 ][2] = data_edge_firstparent_face[ t_e3 ]; // t_e3
//     data_edge_firstparent_face[ t_e3 ] = f3;
//     
//     data_face_nextparents_of_edges[ f4 ][0] = nullindex; // e1
//     data_face_nextparents_of_edges[ f4 ][1] = f5; // e3
//     data_face_nextparents_of_edges[ f4 ][2] = data_edge_firstparent_face[ t_e4 ]; // t_e4
//     data_edge_firstparent_face[ t_e4 ] = f4;
//     
//     data_face_nextparents_of_edges[ f5 ][0] = nullindex; // e2
//     data_face_nextparents_of_edges[ f5 ][1] = nullindex; // e3
//     data_face_nextparents_of_edges[ f5 ][2] = data_edge_firstparent_face[ t_e5 ]; // t_e5
//     data_edge_firstparent_face[ t_e5 ] = f5;
//     
//     data_edge_firstparent_face[ e0 ] = f0;
//     data_edge_firstparent_face[ e1 ] = f0;
//     data_edge_firstparent_face[ e2 ] = f1;
//     data_edge_firstparent_face[ e3 ] = f4;
//     
//     
//     data_edge_nextparents_of_vertices[ e0 ][0] = t_v0;
//     data_edge_nextparents_of_vertices[ e0 ][1] = vn;
//     
//     data_edge_nextparents_of_vertices[ e1 ][0] = t_v1;
//     data_edge_nextparents_of_vertices[ e1 ][1] = vn;
//     
//     data_edge_nextparents_of_vertices[ e2 ][0] = vn;
//     data_edge_nextparents_of_vertices[ e2 ][1] = t_v2;
//     
//     data_edge_nextparents_of_vertices[ e3 ][0] = vn;
//     data_edge_nextparents_of_vertices[ e3 ][1] = t_v2;
//     
//     
//     /* face t_f0: nothing needs to change */
//     
//     /* face t_f1: relink */
//     if( data_face_firstparent_tetrahedron[ t_f1 ] == t ) 
//       data_face_firstparent_tetrahedron[ t_f1 ] = t1;
//     else {
//       int current_tetrahedron = data_face_firstparent_tetrahedron[ t_f1 ];
//       while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] != t 
//              &&
//              data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] != nullindex )
//         current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ];
//       assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] != nullindex );
//       assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] == t );
//       data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f1 ) ] = t1;
//     }
//     
//     /* face t_f2: relink */
//     if( data_face_firstparent_tetrahedron[ t_f2 ] == t ) 
//       data_face_firstparent_tetrahedron[ t_f2 ] = t2;
//     else {
//       int current_tetrahedron = data_face_firstparent_tetrahedron[ t_f2 ];
//       while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] != t 
//              &&
//              data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] != nullindex )
//         current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ];
//       assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] != nullindex );
//       assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] == t );
//       data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f2 ) ] = t2;
//     }
//     
//     /* face t_f3: relink */
//     if( data_face_firstparent_tetrahedron[ t_f3 ] == t ) 
//       data_face_firstparent_tetrahedron[ t_f3 ] = t3;
//     else {
//       int current_tetrahedron = data_face_firstparent_tetrahedron[ t_f3 ];
//       while( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] != t 
//              &&
//              data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] != nullindex )
//         current_tetrahedron = data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ];
//       assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] != nullindex );
//       assert( data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] == t );
//       data_tetrahedron_nextparents_of_faces[ current_tetrahedron ][ indexof_tetrahedron_face( current_tetrahedron, t_f3 ) ] = t3;
//     }
//     
//     
//     /* edge t_e0: nothing needs to change */
//     
//     /* edge t_e1: nothing needs to change */
//     
//     /* edge t_e2: relink */
//     if( data_edge_firstparent_face[ t_e2 ] == t ) 
//       data_edge_firstparent_face[ t_e2 ] = t1; // FIXME
//     else {
//       int current_face = data_edge_firstparent_face[ t_e2 ];
//       while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] != t 
//              &&
//              data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] != nullindex )
//         current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ];
//       assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] != nullindex );
//       assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] == t );
//       data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e2 ) ] = t1; // FIXME
//     }
//         
//     /* edge t_e3: nothing needs to change */
//     
//     /* edge t_e4: relink */
//     if( data_edge_firstparent_face[ t_e4 ] == t ) 
//       data_edge_firstparent_face[ t_e4 ] = t4; // FIXME
//     else {
//       int current_face = data_edge_firstparent_face[ t_e4 ];
//       while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] != t 
//              &&
//              data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] != nullindex )
//         current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ];
//       assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] != nullindex );
//       assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] == t );
//       data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e4 ) ] = t4; // FIXME
//     }
//     
//     /* edge t_e5: relink */
//     if( data_edge_firstparent_face[ t_e5 ] == t ) 
//       data_edge_firstparent_face[ t_e5 ] = t5; // FIXME
//     else {
//       int current_face = data_edge_firstparent_face[ t_e5 ];
//       while( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] != t 
//              &&
//              data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] != nullindex )
//         current_face = data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ];
//       assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] != nullindex );
//       assert( data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] == t );
//       data_face_nextparents_of_edges[ current_face ][ indexof_face_edge( current_face, t_e5 ) ] = t5; // FIXME
//     }
//     
//     /* vertex t_v0: nothing needs to change */
//     
//     /* vertex t_v1: nothing needs to change */
//     
//     /* vertex t_v2: nothing needs to change */
//     
//     /* vertex t_v3: relink */ 
//     if( data_vertex_firstparent_face[ t_v3 ] == t ) 
//       data_vertex_firstparent_face[ t_v3 ] = t1; // FIXME
//     else {
//       int current_face = data_vertex_firstparent_face[ t_v2 ];
//       while( data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] != t 
//              &&
//              data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] != nullindex )
//         current_face = data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ];
//       assert( data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] != nullindex );
//       assert( data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] == t );
//       data_face_nextparents_of_vertices[ current_face ][ indexof_face_vertex( current_face, t_v2 ) ] = t1;
//     }
//     
//     
//     /* face f0: link */
//     data_face_firstparent_tetrahedron[ f0 ] = t0; // FIXME
//     
//     /* face f1: link */
//     data_face_firstparent_tetrahedron[ f1 ] = t0; // FIXME
//     
//     /* face f2: link */
//     data_face_firstparent_tetrahedron[ f2 ] = t0; // FIXME
//     
//     /* face f3: link */
//     data_face_firstparent_tetrahedron[ f3 ] = t0; // FIXME
//     
//     /* face f4: link */
//     data_face_firstparent_tetrahedron[ f4 ] = t0; // FIXME
//     
//     /* face f5: link */
//     data_face_firstparent_tetrahedron[ f5 ] = t0; // FIXME
//     
//     
//     /* edge e0: link */
//     data_edge_firstparent_tetrahedron[ e0 ] = t0; // FIXME
//     data_edge_firstparent_face[ e0 ] = t0; // FIXME
//     
//     /* edge e1: link */
//     data_edge_firstparent_tetrahedron[ e1 ] = t0; // FIXME
//     data_edge_firstparent_face[ e1 ] = t0; // FIXME
//     
//     /* edge e2: link */
//     data_edge_firstparent_tetrahedron[ e2 ] = t0; // FIXME
//     data_edge_firstparent_face[ e2 ] = t1; // FIXME
//     
//     /* edge e3: link */
//     data_edge_firstparent_tetrahedron[ e3 ] = t0; // FIXME
//     data_edge_firstparent_face[ e3 ] = t1; // FIXME
//     
//     
//     /* vertex vn: link */
//     data_vertex_firstparent_tetrahedron[ counter_vertices ] = t0;
//     data_vertex_firstparent_face       [ counter_vertices ] = f0;
//     data_vertex_firstparent_edge       [ counter_vertices ] = e0;
//     
//     
//     /* Counters */
//     counter_tetrahedra = counter_tetrahedra + 3;
//     counter_faces      = counter_faces      + 6;
//     counter_edges      = counter_edges      + 4;
//     counter_vertices   = counter_vertices   + 1;
//     
//     /* DONE */
//     
//     check();
}


void MeshSimplicial3D::midpoint_refinement_global()
{
    check();
    
    int N = counter_tetrahedra;
    
    for( int t = 0; t < N; t++ ) {
      
      midpoint_refinement( t );
      
    }
    
    check();
}











