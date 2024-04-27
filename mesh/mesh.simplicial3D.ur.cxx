
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





void MeshSimplicial3D::uniformrefinement( int levels )
{
  assert( levels >= 0 );
  for( int l = 0; l < levels; l++ )
    uniformrefinement();
}






/*
 *
 *  Split Tetrahedra: 8T 
 *  
 *    00 11 22 33
 *       00 01 02 03
 *       01 11 12 13
 *       02 12 22 23
 *       03 13 23 33
 *       
 *       01 02 03 13
 *       01 02 12 13 
 *       02 03 13 23
 *       02 12 13 23
 *
 *  Split Faces:    4F
 *    
 *    00 11 22
 *       00 01 02
 *       01 11 12
 *       02 12 22
 *       01 02 12
 *    00 11 33
 *       00 01 03
 *       01 11 13
 *       03 13 33
 *       01 03 13
 *    00 22 33
 *       00 02 03
 *       02 22 23
 *       03 23 33
 *       02 03 23
 *    11 22 33
 *       11 12 13
 *       12 22 23
 *       13 23 33
 *       12 13 23
 *       
 *  Completely New Faces: 8T
 *       
 *       01 02 03
 *       01 12 13
 *       02 12 23
 *       03 13 23
 *       
 *       01 02 13
 *       02 03 13
 *       02 12 13
 *       02 13 23
 *       
 * 
 *  Split Edges: 2 E 
 *    
 *    00 11
 *       00 01
 *       01 11
 *    00 22
 *       00 02
 *       02 22
 *    00 33
 *       00 03
 *       03 33
 *    11 22
 *       11 12
 *       12 22
 *    11 33
 *       11 13
 *       13 33
 *    22 33
 *       22 23
 *       23 33
 *    
 *  Completely New Edge: 3F + 1T
 *       
 *       
 *       01 02
 *       01 12
 *       02 12
 *       
 *       01 03
 *       01 13
 *       03 13
 *       
 *       02 03
 *       02 23
 *       03 23
 *       
 *       12 13
 *       12 23
 *       13 23
 *       
 * 
 * 
 *       02 13
 *      
 *
 *
 *
 */




// TODO: Debug the uniform refinement method 

void MeshSimplicial3D::uniformrefinement()
{
    check();
    
    
    int new_counter_tetrahedra = 8 * counter_tetrahedra;
    int new_counter_faces      = 4 * counter_faces + 8 * counter_tetrahedra;
    int new_counter_edges      = 2 * counter_edges + 3  * counter_faces + 1 * counter_tetrahedra;
    int new_counter_vertices   = 1 * counter_vertices + 1 * counter_edges;
    
    
    
    
    
    
    /***********************************/
    /***********************************/
    /******** resize the arrays ********/
    /***********************************/
    /***********************************/
    
    
    
    
    /* tetrahedron -> face */
    
    data_tetrahedron_faces.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    data_face_firstparent_tetrahedron.resize( new_counter_faces, nullindex );
    
    data_tetrahedron_nextparents_of_faces.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
      
    
    /* tetrahedron -> edge */
    
    data_tetrahedron_edges.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    data_edge_firstparent_tetrahedron.resize( new_counter_edges, nullindex );
    
    data_tetrahedron_nextparents_of_edges.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    
    /* tetrahedron -> vertex */
    
    data_tetrahedron_vertices.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    data_vertex_firstparent_tetrahedron.resize( new_counter_vertices, nullindex );
    
    data_tetrahedron_nextparents_of_vertices.resize( new_counter_tetrahedra, { nullindex, nullindex, nullindex, nullindex } );
    
    
    /* face -> edge */
    
    data_face_edges.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    data_edge_firstparent_face.resize( new_counter_edges, nullindex );
    
    data_face_nextparents_of_edges.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    
    /* face -> vertex */
    
    data_face_vertices.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    data_vertex_firstparent_face.resize( new_counter_vertices, nullindex );
    
    data_face_nextparents_of_vertices.resize( new_counter_faces, { nullindex, nullindex, nullindex } );
    
    
    /* edge -> vertex */
    
    data_edge_vertices.resize( new_counter_edges, { nullindex, nullindex } );
    
    data_vertex_firstparent_edge.resize( new_counter_vertices, nullindex );
    
    data_edge_nextparents_of_vertices.resize( new_counter_edges, { nullindex, nullindex } );
    
    
    /* coordinates */
    
    getcoordinates().addcoordinates( counter_edges );
    
    
    
    
    
    
    /* 0. create the new coordinates and fill them up */
    
    for( int e = 0; e < counter_edges; e++ )
    {
      getcoordinates().loadvector( counter_vertices + e, get_edge_midpoint( e ) );
    }
    
    
    
    /***********************************/
    /***********************************/
    /******** set up simplex data ******/
    /***********************************/
    /***********************************/
    
    
    
    /********************************/
    /***   VERTICES AND EDGES    ****/
    /********************************/
    
    /* 1. for each edge created from an old edge, set the vertices */
    
    for( int e = 0; e < counter_edges; e++ )
    {
    
      int vertex_back  = data_edge_vertices[e][0];
      int vertex_front = data_edge_vertices[e][1];
      
      data_edge_vertices[e                ][0] = vertex_back;
      data_edge_vertices[e                ][1] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][0] = counter_vertices + e;
      data_edge_vertices[e + counter_edges][1] = vertex_front;
    
    }
    
    /* 2. for each face, set the vertices of the new edges inside the face */
    
    for( int f = 0; f < counter_faces; f++ )
    {
      
      // 01 02
      data_edge_vertices[ 2 * counter_edges + 0 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][0];
      data_edge_vertices[ 2 * counter_edges + 0 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][1];
      
      // 01 12
      data_edge_vertices[ 2 * counter_edges + 1 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][0];
      data_edge_vertices[ 2 * counter_edges + 1 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][2];
      
      // 02 12
      data_edge_vertices[ 2 * counter_edges + 2 * counter_faces + f ][0] = counter_vertices + data_face_edges[f][1];
      data_edge_vertices[ 2 * counter_edges + 2 * counter_faces + f ][1] = counter_vertices + data_face_edges[f][2];
      
    }
    
    /* 3. for each tetrahedron, set the internal vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      /* 01 02 03 12 13 23 */
      
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int ne = 2 * counter_edges + 3 * counter_faces + t;
      
      data_edge_vertices[ ne ][ 0 ] = v02;
      data_edge_vertices[ ne ][ 1 ] = v13;
      
    }
    
    
    /****************************************/
    /***   VERTICES AND EDGES AND FACES  ****/
    /****************************************/
    
    /*** SET THE VERTICES AND EDGES OF ALL FACES ****/
    
    /* for each new outer face, set the new vertices */
    
    for( int f = 0; f < counter_faces; f++ )
    {
    
      int v00 = data_face_vertices[f][0];
      int v11 = data_face_vertices[f][1];
      int v22 = data_face_vertices[f][2];
      
      int v01 = counter_vertices + data_face_edges[f][0];
      int v02 = counter_vertices + data_face_edges[f][1];
      int v12 = counter_vertices + data_face_edges[f][2];
      
      int f_00_01_12 = 0 * counter_faces + f;
      int f_01_11_12 = 1 * counter_faces + f;
      int f_02_12_22 = 2 * counter_faces + f;
      int f_01_02_12 = 3 * counter_faces + f;
    
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_face_vertices[ f_00_01_12 ][0] = v00;
      data_face_vertices[ f_00_01_12 ][1] = v01;
      data_face_vertices[ f_00_01_12 ][2] = v02;
      
      data_face_vertices[ f_01_11_12 ][0] = v01;
      data_face_vertices[ f_01_11_12 ][1] = v11;
      data_face_vertices[ f_01_11_12 ][2] = v12;
      
      data_face_vertices[ f_02_12_22 ][0] = v02;
      data_face_vertices[ f_02_12_22 ][1] = v12;
      data_face_vertices[ f_02_12_22 ][2] = v22;
      
      data_face_vertices[ f_01_02_12 ][0] = v01;
      data_face_vertices[ f_01_02_12 ][1] = v02;
      data_face_vertices[ f_01_02_12 ][2] = v12;
      
    }
    
    /* for each new interior face, set the new vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int v01 = counter_vertices + data_tetrahedron_edges[t][0];
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v03 = counter_vertices + data_tetrahedron_edges[t][2];
      int v12 = counter_vertices + data_tetrahedron_edges[t][3];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int v23 = counter_vertices + data_tetrahedron_edges[t][5];
      
      assert( counter_vertices <= v01 && v01 < counter_vertices + counter_edges );
      assert( counter_vertices <= v02 && v02 < counter_vertices + counter_edges );
      assert( counter_vertices <= v03 && v03 < counter_vertices + counter_edges );
      assert( counter_vertices <= v12 && v12 < counter_vertices + counter_edges );
      assert( counter_vertices <= v13 && v13 < counter_vertices + counter_edges );
      assert( counter_vertices <= v23 && v23 < counter_vertices + counter_edges );
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      data_face_vertices[ f_01_02_03 ][0] = v01;
      data_face_vertices[ f_01_02_03 ][1] = v02;
      data_face_vertices[ f_01_02_03 ][2] = v03;
      
      data_face_vertices[ f_01_12_13 ][0] = v01;
      data_face_vertices[ f_01_12_13 ][1] = v12;
      data_face_vertices[ f_01_12_13 ][2] = v13;
      
      data_face_vertices[ f_02_12_23 ][0] = v02;
      data_face_vertices[ f_02_12_23 ][1] = v12;
      data_face_vertices[ f_02_12_23 ][2] = v23;
      
      data_face_vertices[ f_03_13_23 ][0] = v03;
      data_face_vertices[ f_03_13_23 ][1] = v13;
      data_face_vertices[ f_03_13_23 ][2] = v23;
      
      data_face_vertices[ f_01_02_13 ][0] = v01;
      data_face_vertices[ f_01_02_13 ][1] = v02;
      data_face_vertices[ f_01_02_13 ][2] = v13;
      
      data_face_vertices[ f_02_03_13 ][0] = v02;
      data_face_vertices[ f_02_03_13 ][1] = v03;
      data_face_vertices[ f_02_03_13 ][2] = v13;
      
      data_face_vertices[ f_02_12_13 ][0] = v02;
      data_face_vertices[ f_02_12_13 ][1] = v12;
      data_face_vertices[ f_02_12_13 ][2] = v13;
      
      data_face_vertices[ f_02_13_23 ][0] = v02;
      data_face_vertices[ f_02_13_23 ][1] = v13;
      data_face_vertices[ f_02_13_23 ][2] = v23;
      
    }
    
    /* for each new outer face, set the new edges */
    
    for( int f = 0; f < counter_faces; f++ )
    {
    
      int e_00_01 = 0 * counter_edges + data_face_edges[f][0];
      int e_00_02 = 0 * counter_edges + data_face_edges[f][1];
      int e_01_11 = 1 * counter_edges + data_face_edges[f][0];
      int e_02_22 = 1 * counter_edges + data_face_edges[f][1];
      int e_11_12 = 0 * counter_edges + data_face_edges[f][2];
      int e_12_22 = 1 * counter_edges + data_face_edges[f][2];
      
      int e_01_02 = 2 * counter_edges + 0 * counter_faces + f;
      int e_01_12 = 2 * counter_edges + 1 * counter_faces + f;
      int e_02_12 = 2 * counter_edges + 2 * counter_faces + f;
      
      // [ 00 01 02 ], [ 01 11 12 ], [ 02 12 22 ], [ 01 02 12 ]
        
      data_face_edges[ 0 * counter_faces + f ][0] = e_00_01;
      data_face_edges[ 0 * counter_faces + f ][1] = e_00_02;
      data_face_edges[ 0 * counter_faces + f ][2] = e_01_02;
      
      data_face_edges[ 1 * counter_faces + f ][0] = e_01_11;
      data_face_edges[ 1 * counter_faces + f ][1] = e_01_12;
      data_face_edges[ 1 * counter_faces + f ][2] = e_11_12;
      
      data_face_edges[ 2 * counter_faces + f ][0] = e_02_12;
      data_face_edges[ 2 * counter_faces + f ][1] = e_02_22;
      data_face_edges[ 2 * counter_faces + f ][2] = e_12_22;
      
      data_face_edges[ 3 * counter_faces + f ][0] = e_01_02;
      data_face_edges[ 3 * counter_faces + f ][1] = e_01_12;
      data_face_edges[ 3 * counter_faces + f ][2] = e_02_12;
      
    }
    
    /* for each new interior face, set the new edges */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      int f_012 = data_tetrahedron_faces[t][0];
      int f_013 = data_tetrahedron_faces[t][1];
      int f_023 = data_tetrahedron_faces[t][2];
      int f_123 = data_tetrahedron_faces[t][3];
      
      int e_01_02 = 2 * counter_edges + 0 * counter_faces + f_012;
      int e_01_12 = 2 * counter_edges + 1 * counter_faces + f_012;
      int e_02_12 = 2 * counter_edges + 2 * counter_faces + f_012;
      
      int e_01_03 = 2 * counter_edges + 0 * counter_faces + f_013;
      int e_01_13 = 2 * counter_edges + 1 * counter_faces + f_013;
      int e_03_13 = 2 * counter_edges + 2 * counter_faces + f_013;
      
      int e_02_03 = 2 * counter_edges + 0 * counter_faces + f_023;
      int e_02_23 = 2 * counter_edges + 1 * counter_faces + f_023;
      int e_03_23 = 2 * counter_edges + 2 * counter_faces + f_023;
      
      int e_12_13 = 2 * counter_edges + 0 * counter_faces + f_123;
      int e_12_23 = 2 * counter_edges + 1 * counter_faces + f_123;
      int e_13_23 = 2 * counter_edges + 2 * counter_faces + f_123;
      
      int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
      
      data_face_edges[ f_01_02_03 ][0] = e_01_02;
      data_face_edges[ f_01_02_03 ][1] = e_01_03;
      data_face_edges[ f_01_02_03 ][2] = e_02_03;
      
      data_face_edges[ f_01_12_13 ][0] = e_01_12;
      data_face_edges[ f_01_12_13 ][1] = e_01_13;
      data_face_edges[ f_01_12_13 ][2] = e_12_13;
      
      data_face_edges[ f_02_12_23 ][0] = e_02_12;
      data_face_edges[ f_02_12_23 ][1] = e_02_23;
      data_face_edges[ f_02_12_23 ][2] = e_12_23;
      
      data_face_edges[ f_03_13_23 ][0] = e_03_13;
      data_face_edges[ f_03_13_23 ][1] = e_03_23;
      data_face_edges[ f_03_13_23 ][2] = e_13_23;
      
      data_face_edges[ f_01_02_13 ][0] = e_01_02;
      data_face_edges[ f_01_02_13 ][1] = e_01_13;
      data_face_edges[ f_01_02_13 ][2] = e_02_13;
      
      data_face_edges[ f_02_03_13 ][0] = e_02_03;
      data_face_edges[ f_02_03_13 ][1] = e_02_13;
      data_face_edges[ f_02_03_13 ][2] = e_03_13;
      
      data_face_edges[ f_02_12_13 ][0] = e_02_12;
      data_face_edges[ f_02_12_13 ][1] = e_02_13;
      data_face_edges[ f_02_12_13 ][2] = e_12_13;
      
      data_face_edges[ f_02_13_23 ][0] = e_02_13;
      data_face_edges[ f_02_13_23 ][1] = e_02_23;
      data_face_edges[ f_02_13_23 ][2] = e_13_23;
      
    }
    
    
    
    /**************************/
    /***   ALL DIMENSIONS  ****/
    /**************************/
    
    /*
    *       00 01 02 03
    *       01 11 12 13
    *       02 12 22 23
    *       03 13 23 33
    *       
    *       01 02 03 13
    *       01 02 12 13 
    *       02 03 13 23
    *       02 12 13 23
    */
    
    /* for each new tetrahedron, set the new vertices */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      
      int v00 = data_tetrahedron_vertices[t][0];
      int v11 = data_tetrahedron_vertices[t][1];
      int v22 = data_tetrahedron_vertices[t][2];
      int v33 = data_tetrahedron_vertices[t][3];
      
      int v01 = counter_vertices + data_tetrahedron_edges[t][0];
      int v02 = counter_vertices + data_tetrahedron_edges[t][1];
      int v03 = counter_vertices + data_tetrahedron_edges[t][2];
      int v12 = counter_vertices + data_tetrahedron_edges[t][3];
      int v13 = counter_vertices + data_tetrahedron_edges[t][4];
      int v23 = counter_vertices + data_tetrahedron_edges[t][5];
      
      //       00 01 02 03
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][0] = v00;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][1] = v01;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][2] = v02;
      data_tetrahedron_vertices[ 0 * counter_tetrahedra + t ][3] = v03;
      
      //       01 11 12 13
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][1] = v11;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][2] = v12;
      data_tetrahedron_vertices[ 1 * counter_tetrahedra + t ][3] = v13;
      
      //       02 12 22 23
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][1] = v12;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][2] = v22;
      data_tetrahedron_vertices[ 2 * counter_tetrahedra + t ][3] = v23;
      
      //       03 13 23 33
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][0] = v03;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][1] = v13;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][2] = v23;
      data_tetrahedron_vertices[ 3 * counter_tetrahedra + t ][3] = v33;
      
      
      //       01 02 03 13
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][1] = v02;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][2] = v03;
      data_tetrahedron_vertices[ 4 * counter_tetrahedra + t ][3] = v13;
      
      //       01 02 12 13 
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][0] = v01;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][1] = v02;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][2] = v12;
      data_tetrahedron_vertices[ 5 * counter_tetrahedra + t ][3] = v13;
      
      //       02 03 13 23
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][1] = v03;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][2] = v13;
      data_tetrahedron_vertices[ 6 * counter_tetrahedra + t ][3] = v23;
      
      //       02 12 13 23
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][0] = v02;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][1] = v12;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][2] = v13;
      data_tetrahedron_vertices[ 7 * counter_tetrahedra + t ][3] = v23;
      
    }
    
    /* for each new tetrahedron, set the new edges */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
        // 00 11 22 33 
        
        // From old edges ...
        int e_00_01 = data_tetrahedron_edges[t][0] + 0 * counter_edges;
        int e_01_11 = data_tetrahedron_edges[t][0] + 1 * counter_edges;
        
        int e_00_02 = data_tetrahedron_edges[t][1] + 0 * counter_edges;
        int e_02_22 = data_tetrahedron_edges[t][1] + 1 * counter_edges;
        
        int e_00_03 = data_tetrahedron_edges[t][2] + 0 * counter_edges;
        int e_03_33 = data_tetrahedron_edges[t][2] + 1 * counter_edges;
        
        int e_11_12 = data_tetrahedron_edges[t][3] + 0 * counter_edges;
        int e_12_22 = data_tetrahedron_edges[t][3] + 1 * counter_edges;
        
        int e_11_13 = data_tetrahedron_edges[t][4] + 0 * counter_edges;
        int e_13_33 = data_tetrahedron_edges[t][4] + 1 * counter_edges;
        
        int e_22_23 = data_tetrahedron_edges[t][5] + 0 * counter_edges;
        int e_23_33 = data_tetrahedron_edges[t][5] + 1 * counter_edges;
        
        
        // ... from old faces ... 
        
        // 00 11 22  
        int e_01_02 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][0];
        int e_01_12 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][0];
        int e_02_12 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][0];
        
        // 00 11 33 
        int e_01_03 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][1];
        int e_01_13 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][1];
        int e_03_13 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][1];
        
        // 00 22 33 
        int e_02_03 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][2];
        int e_02_23 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][2];
        int e_03_23 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][2];
        
        // 11 22 33 
        int e_12_13 = 2 * counter_edges + 0 * counter_faces + data_tetrahedron_faces[t][3];
        int e_12_23 = 2 * counter_edges + 1 * counter_faces + data_tetrahedron_faces[t][3];
        int e_13_23 = 2 * counter_edges + 2 * counter_faces + data_tetrahedron_faces[t][3];
        
        // ... and the single new internal one. 
        
        int e_02_13 = 2 * counter_edges + 3 * counter_faces + t;
        
        
        // fill in to the tetrahedra
        
        //       00 01 02 03
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][0] = e_00_01;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][1] = e_00_02;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][2] = e_00_03;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][3] = e_01_02;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][4] = e_01_03;
        data_tetrahedron_edges[ 0 * counter_tetrahedra + t ][5] = e_02_03;
        
        //       01 11 12 13
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][0] = e_01_11;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][1] = e_01_12;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][2] = e_01_13;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][3] = e_11_12;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][4] = e_11_13;
        data_tetrahedron_edges[ 1 * counter_tetrahedra + t ][5] = e_12_13;
        
        //       02 12 22 23
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][0] = e_02_12;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][1] = e_02_22;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][2] = e_02_23;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][3] = e_12_22;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][4] = e_12_23;
        data_tetrahedron_edges[ 2 * counter_tetrahedra + t ][5] = e_22_23;
        
        //       03 13 23 33
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][0] = e_03_13;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][1] = e_03_23;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][2] = e_03_33;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][3] = e_13_23;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][4] = e_13_33;
        data_tetrahedron_edges[ 3 * counter_tetrahedra + t ][5] = e_23_33;
        
        
        //       01 02 03 13
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][0] = e_01_02;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][1] = e_01_03;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][2] = e_01_13;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][3] = e_02_03;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][4] = e_02_13;
        data_tetrahedron_edges[ 4 * counter_tetrahedra + t ][5] = e_03_13;
        
        //       01 02 12 13 
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][0] = e_01_02;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][1] = e_01_12;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][2] = e_01_13;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][3] = e_02_12;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][4] = e_02_13;
        data_tetrahedron_edges[ 5 * counter_tetrahedra + t ][5] = e_12_13;
        
        //       02 03 13 23
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][0] = e_02_03;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][1] = e_02_13;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][2] = e_02_23;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][3] = e_03_13;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][4] = e_03_23;
        data_tetrahedron_edges[ 6 * counter_tetrahedra + t ][5] = e_13_23;
        
        //       02 12 13 23
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][0] = e_02_12;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][1] = e_02_13;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][2] = e_02_23;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][3] = e_12_13;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][4] = e_12_23;
        data_tetrahedron_edges[ 7 * counter_tetrahedra + t ][5] = e_13_23;
        
    }
    
    /* for each new tetrahedron, set the new faces */
    
    for( int t = 0; t < counter_tetrahedra; t++ )
    {
      // 00 11 22 33 
      
      // From old faces ... 
      
      // 00 11 22 
      int f_00_01_02 = data_tetrahedron_faces[t][0] + 0 * counter_faces;
      int f_01_11_12 = data_tetrahedron_faces[t][0] + 1 * counter_faces;
      int f_02_12_22 = data_tetrahedron_faces[t][0] + 2 * counter_faces;
      int f_01_02_12 = data_tetrahedron_faces[t][0] + 3 * counter_faces;
      
      // 00 11 33 
      int f_00_01_03 = data_tetrahedron_faces[t][1] + 0 * counter_faces;
      int f_01_11_13 = data_tetrahedron_faces[t][1] + 1 * counter_faces;
      int f_03_13_33 = data_tetrahedron_faces[t][1] + 2 * counter_faces;
      int f_01_03_13 = data_tetrahedron_faces[t][1] + 3 * counter_faces;
      
      // 00 22 33 
      int f_00_02_03 = data_tetrahedron_faces[t][2] + 0 * counter_faces;
      int f_02_22_23 = data_tetrahedron_faces[t][2] + 1 * counter_faces;
      int f_03_23_33 = data_tetrahedron_faces[t][2] + 2 * counter_faces;
      int f_02_03_23 = data_tetrahedron_faces[t][2] + 3 * counter_faces;
      
      // 11 22 33 
      int f_11_12_13 = data_tetrahedron_faces[t][3] + 0 * counter_faces;
      int f_12_22_23 = data_tetrahedron_faces[t][3] + 1 * counter_faces;
      int f_13_23_33 = data_tetrahedron_faces[t][3] + 2 * counter_faces;
      int f_12_13_23 = data_tetrahedron_faces[t][3] + 3 * counter_faces;
      
      // ... new internal faces. 
      
      // 01 02 03
      // 01 12 13
      // 02 12 23
      // 03 13 23
      int f_01_02_03 = 4 * counter_faces + 0 * counter_tetrahedra + t;
      int f_01_12_13 = 4 * counter_faces + 1 * counter_tetrahedra + t;
      int f_02_12_23 = 4 * counter_faces + 2 * counter_tetrahedra + t;
      int f_03_13_23 = 4 * counter_faces + 3 * counter_tetrahedra + t;
      
      // 01 02 13
      // 02 03 13
      // 02 12 13
      // 02 13 23
      int f_01_02_13 = 4 * counter_faces + 4 * counter_tetrahedra + t;
      int f_02_03_13 = 4 * counter_faces + 5 * counter_tetrahedra + t;
      int f_02_12_13 = 4 * counter_faces + 6 * counter_tetrahedra + t;
      int f_02_13_23 = 4 * counter_faces + 7 * counter_tetrahedra + t;
      
      
      
      
      //       00 01 02 03
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][0] = f_00_01_02;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][1] = f_00_01_03;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][2] = f_00_02_03;
      data_tetrahedron_faces[ 0 * counter_tetrahedra + t ][3] = f_01_02_03;
      
      //       01 11 12 13
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][0] = f_01_11_12;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][1] = f_01_11_13;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][2] = f_01_12_13;
      data_tetrahedron_faces[ 1 * counter_tetrahedra + t ][3] = f_11_12_13;
      
      //       02 12 22 23
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][0] = f_02_12_22;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][1] = f_02_12_23;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][2] = f_02_22_23;
      data_tetrahedron_faces[ 2 * counter_tetrahedra + t ][3] = f_12_22_23;
      
      //       03 13 23 33
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][0] = f_03_13_23;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][1] = f_03_13_33;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][2] = f_03_23_33;
      data_tetrahedron_faces[ 3 * counter_tetrahedra + t ][3] = f_13_23_33;
      
      
      //       01 02 03 13
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][0] = f_01_02_03;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][1] = f_01_02_13;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][2] = f_01_03_13;
      data_tetrahedron_faces[ 4 * counter_tetrahedra + t ][3] = f_02_03_13;
      
      //       01 02 12 13 
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][0] = f_01_02_12;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][1] = f_01_02_13;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][2] = f_01_12_13;
      data_tetrahedron_faces[ 5 * counter_tetrahedra + t ][3] = f_02_12_13;
      
      //       02 03 13 23
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][0] = f_02_03_13;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][1] = f_02_03_23;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][2] = f_02_13_23;
      data_tetrahedron_faces[ 6 * counter_tetrahedra + t ][3] = f_03_13_23;
      
      //       02 12 13 23
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][0] = f_02_12_13;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][1] = f_02_12_23;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][2] = f_02_13_23;
      data_tetrahedron_faces[ 7 * counter_tetrahedra + t ][3] = f_12_13_23;
      
    }
    
    
    
    
    
    
    /***********************************/
    /***********************************/
    /**** set up list of parents *******/
    /***********************************/
    /***********************************/
    
    
    std::fill( data_vertex_firstparent_edge.begin(), data_vertex_firstparent_edge.end(), nullindex );
    
    std::fill( data_edge_nextparents_of_vertices.begin(), data_edge_nextparents_of_vertices.end(), std::array<int,2>{ nullindex, nullindex } );
    
    for( int e = 0; e < new_counter_edges; e++ )
    for( int vi = 0; vi < 2; vi++ )
    {
        
        int v = data_edge_vertices[ e ][ vi ];
        int parent = data_vertex_firstparent_edge[ v ];
        data_vertex_firstparent_edge[ v ] = e;
        data_edge_nextparents_of_vertices[ e ][ vi ] = parent;
        
    }
    
    
    std::fill( data_vertex_firstparent_face.begin(), data_vertex_firstparent_face.end(), nullindex );
    
    std::fill( data_face_nextparents_of_vertices.begin(), data_face_nextparents_of_vertices.end(), std::array<int,3>{ nullindex, nullindex, nullindex } );

    for( int f = 0; f < new_counter_faces; f++ )
    for( int vi = 0; vi < 3; vi++ )
    {
        
        int v = data_face_vertices[ f ][ vi ];
        int parent = data_vertex_firstparent_face[ v ];
        data_vertex_firstparent_face[ v ] = f;
        data_face_nextparents_of_vertices[ f ][ vi ] = parent;
        
    }
    
    
    std::fill( data_edge_firstparent_face.begin(), data_edge_firstparent_face.end(), nullindex );
    
    std::fill( data_face_nextparents_of_edges.begin(), data_face_nextparents_of_edges.end(), std::array<int,3>{ nullindex, nullindex, nullindex } );
    
    for( int f = 0; f < new_counter_faces; f++ )
    for( int ei = 0; ei < 3; ei++ )
    {
        
        int e = data_face_edges[ f ][ ei ];
        int parent = data_edge_firstparent_face[ e ];
        data_edge_firstparent_face[ e ] = f;
        data_face_nextparents_of_edges[ f ][ ei ] = parent;
        
    }
    
    
    std::fill( data_vertex_firstparent_tetrahedron.begin(), data_vertex_firstparent_tetrahedron.end(), nullindex );
     
    std::fill( data_tetrahedron_nextparents_of_vertices.begin(), data_tetrahedron_nextparents_of_vertices.end(), std::array<int,4>{ nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < new_counter_tetrahedra; t++ )
    for( int vi = 0; vi < 4; vi++ )
    {
        
        int v = data_tetrahedron_vertices[ t ][ vi ];
        int parent = data_vertex_firstparent_tetrahedron[ v ];
        data_vertex_firstparent_tetrahedron[ v ] = t;
        data_tetrahedron_nextparents_of_vertices[ t ][ vi ] = parent;
        
    }
    
    
    std::fill( data_edge_firstparent_tetrahedron.begin(), data_edge_firstparent_tetrahedron.end(), nullindex );
    
    std::fill( data_tetrahedron_nextparents_of_edges.begin(), data_tetrahedron_nextparents_of_edges.end(), std::array<int,6>{ nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < new_counter_tetrahedra; t++ )
    for( int ei = 0; ei < 6; ei++ )
    {
        
        int e = data_tetrahedron_edges[ t ][ ei ];
        int parent = data_edge_firstparent_tetrahedron[ e ];
        data_edge_firstparent_tetrahedron[ e ] = t;
        data_tetrahedron_nextparents_of_edges[ t ][ ei ] = parent;
        
    }

    
    std::fill( data_face_firstparent_tetrahedron.begin(), data_face_firstparent_tetrahedron.end(), nullindex );
    
    std::fill( data_tetrahedron_nextparents_of_faces.begin(), data_tetrahedron_nextparents_of_faces.end(), std::array<int,4>{ nullindex, nullindex, nullindex, nullindex } );
    
    for( int t = 0; t < new_counter_tetrahedra; t++ )
    for( int fi = 0; fi < 4; fi++ )
    {
        
        int f = data_tetrahedron_faces[ t ][ fi ];
        int parent = data_face_firstparent_tetrahedron[ f ];
        data_face_firstparent_tetrahedron[ f ] = t;
        data_tetrahedron_nextparents_of_faces[ t ][ fi ] = parent;
        
    }
    
    
    
    flags_tetrahedra.resize( new_counter_tetrahedra, SimplexFlag::SimplexFlagInvalid );
    flags_faces.resize     ( new_counter_faces,      SimplexFlag::SimplexFlagInvalid );
    flags_edges.resize     ( new_counter_edges,      SimplexFlag::SimplexFlagInvalid );
    flags_vertices.resize  ( new_counter_vertices,   SimplexFlag::SimplexFlagInvalid );
    
    // flags of vertices/edges created from edges vertices
    for( int e = 0; e < counter_edges; e++ ) { 
        
        auto flag_of_edge = flags_edges[ e ];
        
        flags_edges[ counter_edges + e ] = flag_of_edge;
        
        flags_vertices[ counter_vertices + e ] = flag_of_edge;
    }
    
    // flags of edges/faces created from faces 
    for( int f = 0; f < counter_faces; f++ ) {
        
        auto flag_of_face = flags_faces[ f ];
        
        flags_faces[ 1 * counter_faces + f ] = flag_of_face;
        flags_faces[ 2 * counter_faces + f ] = flag_of_face;
        flags_faces[ 3 * counter_faces + f ] = flag_of_face;
        
        flags_edges[ 2 * counter_edges + 0 * counter_faces + f ] = flag_of_face;
        flags_edges[ 2 * counter_edges + 1 * counter_faces + f ] = flag_of_face;
        flags_edges[ 2 * counter_edges + 2 * counter_faces + f ] = flag_of_face;
    }

    // flags of edges/faces/tetrahedra created from tetrahedra
    for( int t = 0; t < counter_tetrahedra; t++ ) {
        
        auto flag_of_tet = flags_tetrahedra[ t ];
        
        flags_tetrahedra[ 1 * counter_tetrahedra + t ] = flag_of_tet;
        flags_tetrahedra[ 2 * counter_tetrahedra + t ] = flag_of_tet;
        flags_tetrahedra[ 3 * counter_tetrahedra + t ] = flag_of_tet;
        flags_tetrahedra[ 4 * counter_tetrahedra + t ] = flag_of_tet;
        flags_tetrahedra[ 5 * counter_tetrahedra + t ] = flag_of_tet;
        flags_tetrahedra[ 6 * counter_tetrahedra + t ] = flag_of_tet;
        flags_tetrahedra[ 7 * counter_tetrahedra + t ] = flag_of_tet;
        
        flags_faces[ 4 * counter_faces + 0 * counter_tetrahedra + t ] = flag_of_tet;
        flags_faces[ 4 * counter_faces + 1 * counter_tetrahedra + t ] = flag_of_tet;
        flags_faces[ 4 * counter_faces + 2 * counter_tetrahedra + t ] = flag_of_tet;
        flags_faces[ 4 * counter_faces + 3 * counter_tetrahedra + t ] = flag_of_tet;
        flags_faces[ 4 * counter_faces + 4 * counter_tetrahedra + t ] = flag_of_tet;
        flags_faces[ 4 * counter_faces + 5 * counter_tetrahedra + t ] = flag_of_tet;
        flags_faces[ 4 * counter_faces + 6 * counter_tetrahedra + t ] = flag_of_tet;
        flags_faces[ 4 * counter_faces + 7 * counter_tetrahedra + t ] = flag_of_tet;
        
        flags_edges[ 2 * counter_edges + 3 * counter_faces + t ] = flag_of_tet;
    }

    
    
    
    
    /* update the counters */
    
    counter_vertices   = new_counter_vertices;
    counter_edges      = new_counter_edges;
    counter_faces      = new_counter_faces;
    counter_tetrahedra = new_counter_tetrahedra;
    
    /* DONE */
    
    check();
}








