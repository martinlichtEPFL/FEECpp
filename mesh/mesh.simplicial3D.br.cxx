
#include <string>
#include <vector>
#include <stack> // TODO: change to something else such as list 
#include <list>
#include <map>
#include <set>
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















/*
 * * * * * VERTEX BISECTION
 */
    
void MeshSimplicial3D::bisect_edge( int e )
{
    assert( 0 <= e && e < counter_edges );
    // check();
    
    
    
    
    /*******************/
    /* DATA COLLECTION */
    /*******************/
    
    // gather basic data about the bisection, in particular the number of simplices involved.
    
    int e_back_vertex  = data_edge_vertices[ e ][ 0 ];
    int e_front_vertex = data_edge_vertices[ e ][ 1 ];
    
    std::vector<int> old_faces      = get_face_parents_of_edge( e );
    std::vector<int> old_tetrahedra = get_tetrahedron_parents_of_edge( e );
    
    std::vector<int> localindex_of_face_refinementedge( old_faces.size() );
    std::vector<int> localindex_of_tetrahedron_refinementedge( old_tetrahedra.size() );
    
    
    FloatVector midcoordinate = get_edge_midpoint( e );
    
    
    
    
    /*******************************************/
    /*****                                ******/
    /*****   Assemble sets of simplices   ******/
    /*****   in the micropatch            ******/
    /*****                                ******/
    /*******************************************/
    
    std::set<int> micropatch_tetrahedra;
    std::set<int> micropatch_faces;
    std::set<int> micropatch_edges;
    std::set<int> micropatch_vertices;
    
    for( int t : old_tetrahedra ){
        
        micropatch_tetrahedra.insert( t );
        
        for( int f : data_tetrahedron_faces[t]  )
            micropatch_faces.insert( f );
        
        for( int e : data_tetrahedron_edges[t]  )
            micropatch_edges.insert( e );
        
        for( int v : data_tetrahedron_vertices[t]  )
            micropatch_vertices.insert( v );
                
    }
    
    
    /*******************************************/
    /*****                                ******/
    /*****   Assemble sets of simplices   ******/
    /*****   in the MACROpatch            ******/
    /*****                                ******/
    /*******************************************/
    
    std::set<int> macropatch_tetrahedra;
    std::set<int> macropatch_faces;
    std::set<int> macropatch_edges;
    std::set<int> macropatch_vertices;
    
    for( int v : micropatch_vertices ){
        
        for(
            int t = get_vertex_firstparent_tetrahedron( v );
            t != nullindex; 
            t = get_vertex_nextparent_tetrahedron( v ,t )
        ) {
            
            macropatch_tetrahedra.insert( t );
            
            for( int fi = 0; fi < 4; fi++ )
                macropatch_faces.insert( data_tetrahedron_faces[t][fi] );
            
            for( int ei = 0; ei < 6; ei++ )
                macropatch_edges.insert( data_tetrahedron_edges[t][ei] );
            
            for( int vi = 0; vi < 4; vi++ )
                macropatch_vertices.insert( data_tetrahedron_vertices[t][vi] );
            
        }
        
    }
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /***********************/
    /*                     */
    /*   ALLOCATE MEMORY   */
    /*                     */
    /***********************/
    
    // allocate additional memory for the new simplex data after bisection 
    
    data_tetrahedron_nextparents_of_faces.resize( counter_tetrahedra + old_tetrahedra.size(),               { nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_faces.resize               ( counter_tetrahedra + old_tetrahedra.size(),               { nullindex, nullindex, nullindex, nullindex } );
    data_face_firstparent_tetrahedron.resize    ( counter_faces + old_faces.size() + old_tetrahedra.size(),                                      nullindex );
    
    data_tetrahedron_nextparents_of_edges.resize( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_edges.resize               ( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex, nullindex, nullindex } );
    data_edge_firstparent_tetrahedron.resize    ( counter_edges + old_faces.size() + 1,                                                                  nullindex );
    
    data_tetrahedron_nextparents_of_vertices.resize( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex } );
    data_tetrahedron_vertices.resize               ( counter_tetrahedra + old_tetrahedra.size(), { nullindex, nullindex, nullindex, nullindex } );
    data_vertex_firstparent_tetrahedron.resize     ( counter_vertices + 1,                                                            nullindex );
    
    data_face_nextparents_of_edges.resize( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_face_edges.resize               ( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_edge_firstparent_face.resize    ( counter_edges     + old_faces.size() + 1,                                           nullindex );
    
    data_face_nextparents_of_vertices.resize( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_face_vertices.resize               ( counter_faces + old_faces.size() + old_tetrahedra.size(), { nullindex, nullindex, nullindex } );
    data_vertex_firstparent_face.resize     ( counter_vertices + 1,                                                               nullindex );
    
    data_edge_nextparents_of_vertices.resize( counter_edges    + old_faces.size() + 1, { nullindex, nullindex } );
    data_edge_vertices.resize               ( counter_edges    + old_faces.size() + 1, { nullindex, nullindex } );
    data_vertex_firstparent_edge.resize     ( counter_vertices + 1,                                   nullindex );
    
    
    
    /* flags */
    {
        flags_tetrahedra.resize( counter_tetrahedra + old_tetrahedra.size(),               SimplexFlag::SimplexFlagInvalid );
        flags_faces.resize     ( counter_faces + old_faces.size() + old_tetrahedra.size(), SimplexFlag::SimplexFlagInvalid );
        flags_edges.resize     ( counter_edges + 1 + old_faces.size(),                     SimplexFlag::SimplexFlagInvalid );
        flags_vertices.resize  ( counter_vertices + 1,                                     SimplexFlag::SimplexFlagInvalid );
        
        SimplexFlag flag_e = flags_edges[e];
        flags_vertices[ counter_vertices ] = flag_e;
        flags_edges   [ counter_edges    ] = flag_e;
        
        for( int of = 0; of < old_faces.size(); of++ ){
            SimplexFlag flag_of = flags_faces[ old_faces[of] ];
            flags_faces[ counter_faces + of     ] = flag_of;
            flags_edges[ counter_edges + 1 + of ] = flag_of;
        }
        
        for( int ot = 0; ot < old_tetrahedra.size(); ot++ ){
            SimplexFlag flag_ot = flags_tetrahedra[ old_tetrahedra[ot] ];
            flags_tetrahedra[ counter_tetrahedra + ot               ] = flag_ot;
            flags_faces     [ counter_faces + old_faces.size() + ot ] = flag_ot;
        }
    }
    
    
    
    
    
    /*********************/
    /*                   */
    /*   GEOMETRY DATA   */
    /*                   */
    /*********************/
    
    getcoordinates().append( midcoordinate );
    
    
    
    
    /****************************/
    /*                          */
    /*   FILL IN SIMPLEX DATA   */
    /*                          */
    /****************************/
    
    
    
    /* insert the new vertex into the micropatch and macropatch */
    
    micropatch_vertices.insert( counter_vertices );
    macropatch_vertices.insert( counter_vertices );
    
    
    
    
    /* vertices of the children of the bisected edge */
    
    data_edge_vertices[ e             ][ 0 ] = e_back_vertex;
    data_edge_vertices[ e             ][ 1 ] = counter_vertices;
      
    data_edge_vertices[ counter_edges ][ 0 ] = counter_vertices;
    data_edge_vertices[ counter_edges ][ 1 ] = e_front_vertex;
    
    
    /* insert the new edge into the micropatch and macropatch */
    
    micropatch_edges.insert( counter_edges );
    macropatch_edges.insert( counter_edges );
    
    
    
    
    
    for( int of = 0; of < old_faces.size(); of++ ) {
      
      // gather the previous data from the current old face   
      
      int f_old = old_faces[ of ];
      int f_new = counter_faces + of;
      
      int e_0_1 = data_face_edges[ f_old ][ 0 ];
      int e_0_2 = data_face_edges[ f_old ][ 1 ];
      int e_1_2 = data_face_edges[ f_old ][ 2 ];
      
      int v_0 = data_face_vertices[ f_old ][ 0 ];
      int v_1 = data_face_vertices[ f_old ][ 1 ];
      int v_2 = data_face_vertices[ f_old ][ 2 ];
      
      
      localindex_of_face_refinementedge[ of ] = ( e_0_1 == e ) ? 0 : ( ( e_0_2 == e ) ? 1 : 2 );
      assert( data_face_edges[ f_old ][ localindex_of_face_refinementedge[ of ] ] == e );
      
      
      /* insert the new face and the new edge into the micropatch and macropatch */
      
      micropatch_faces.insert( f_new ); 
      macropatch_faces.insert( f_new ); 
      
      micropatch_edges.insert( counter_edges + 1 + of );
      macropatch_edges.insert( counter_edges + 1 + of );
      
      
      if( localindex_of_face_refinementedge[ of ] == 0 ) { // 0 1 
        
        assert( v_0 == e_back_vertex && v_1 == e_front_vertex && e == e_0_1 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = v_0;
        data_face_vertices[ f_old ][1] = counter_vertices;
        data_face_vertices[ f_old ][2] = v_2;
        
        data_face_vertices[ f_new ][0] = counter_vertices;
        data_face_vertices[ f_new ][1] = v_1;
        data_face_vertices[ f_new ][2] = v_2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = e;
        data_face_edges[ f_old ][1] = e_0_2;
        data_face_edges[ f_old ][2] = counter_edges + 1 + of;
        
        data_face_edges[ f_new ][0] = counter_edges;
        data_face_edges[ f_new ][1] = counter_edges + 1 + of;
        data_face_edges[ f_new ][2] = e_1_2;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = counter_vertices;
        data_edge_vertices[ counter_edges + 1 + of ][1] = v_2;
        
      } else if( localindex_of_face_refinementedge[ of ] == 1 ) { // 0 2 
        
        assert( v_0 == e_back_vertex && v_2 == e_front_vertex && e == e_0_2 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = v_0;
        data_face_vertices[ f_old ][1] = v_1;
        data_face_vertices[ f_old ][2] = counter_vertices;
        
        data_face_vertices[ f_new ][0] = v_1;
        data_face_vertices[ f_new ][1] = counter_vertices;
        data_face_vertices[ f_new ][2] = v_2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = e_0_1;
        data_face_edges[ f_old ][1] = e;
        data_face_edges[ f_old ][2] = counter_edges + 1 + of;
        
        data_face_edges[ f_new ][0] = counter_edges + 1 + of;
        data_face_edges[ f_new ][1] = e_1_2;
        data_face_edges[ f_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = v_1;
        data_edge_vertices[ counter_edges + 1 + of ][1] = counter_vertices;
        
      } else if( localindex_of_face_refinementedge[ of ] == 2 ) { // 1 2 
        
        assert( v_1 == e_back_vertex && v_2 == e_front_vertex && e == e_1_2 );
        
        /* face vertices */
        data_face_vertices[ f_old ][0] = v_0;
        data_face_vertices[ f_old ][1] = v_1;
        data_face_vertices[ f_old ][2] = counter_vertices;
        
        data_face_vertices[ f_new ][0] = v_0;
        data_face_vertices[ f_new ][1] = counter_vertices;
        data_face_vertices[ f_new ][2] = v_2;
        
        /* face edges */
        data_face_edges[ f_old ][0] = e_0_1;
        data_face_edges[ f_old ][1] = counter_edges + 1 + of;
        data_face_edges[ f_old ][2] = e;
        
        data_face_edges[ f_new ][0] = counter_edges + 1 + of;
        data_face_edges[ f_new ][1] = e_0_2;
        data_face_edges[ f_new ][2] = counter_edges;
        
        /* bisection edge vertices */
        data_edge_vertices[ counter_edges + 1 + of ][0] = v_0;
        data_edge_vertices[ counter_edges + 1 + of ][1] = counter_vertices;
        
      } else {
        
        unreachable();
        
      } 
      
    }
    
    
    
    for( int ot = 0; ot < old_tetrahedra.size(); ot++ ) {
      
      int t_old = old_tetrahedra[ ot ];
      int t_new = counter_tetrahedra + ot;
      
      int f_new = counter_faces + old_faces.size() + ot;
      
      int f_012 = data_tetrahedron_faces[ t_old ][ 0 ];
      int f_013 = data_tetrahedron_faces[ t_old ][ 1 ];
      int f_023 = data_tetrahedron_faces[ t_old ][ 2 ];
      int f_123 = data_tetrahedron_faces[ t_old ][ 3 ];
      
      int e_0_1 = data_tetrahedron_edges[ t_old ][ 0 ];
      int e_0_2 = data_tetrahedron_edges[ t_old ][ 1 ];
      int e_0_3 = data_tetrahedron_edges[ t_old ][ 2 ];
      int e_1_2 = data_tetrahedron_edges[ t_old ][ 3 ];
      int e_1_3 = data_tetrahedron_edges[ t_old ][ 4 ];
      int e_2_3 = data_tetrahedron_edges[ t_old ][ 5 ];
      
      int v_0 = data_tetrahedron_vertices[ t_old ][ 0 ];
      int v_1 = data_tetrahedron_vertices[ t_old ][ 1 ];
      int v_2 = data_tetrahedron_vertices[ t_old ][ 2 ];
      int v_3 = data_tetrahedron_vertices[ t_old ][ 3 ];
      
      int v_n = counter_vertices; 
      
      localindex_of_tetrahedron_refinementedge[ ot ] = nullindex;
      if( e_0_1 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 0;
      if( e_0_2 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 1;
      if( e_0_3 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 2;
      if( e_1_2 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 3;
      if( e_1_3 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 4;
      if( e_2_3 == e ) localindex_of_tetrahedron_refinementedge[ ot ] = 5;
      assert( localindex_of_tetrahedron_refinementedge[ ot ] != nullindex );
      assert( data_tetrahedron_edges[ t_old ][ localindex_of_tetrahedron_refinementedge[ ot ] ] == e );
      
      
      /* insert the new face and the new edge into the micropatch and macropatch */
      
      micropatch_tetrahedra.insert( t_new ); 
      macropatch_tetrahedra.insert( t_new ); 
      
      micropatch_faces.insert( f_new ); 
      macropatch_faces.insert( f_new ); 
      
      
      
      
      
      
      if(        localindex_of_tetrahedron_refinementedge[ ot ] == 0 ) { // 0 1 
        
        assert( v_0 == e_back_vertex && v_1 == e_front_vertex && e == e_0_1 );
        
        /* prepare stuff */ 
        int index_f_012 = std::find( old_faces.begin(), old_faces.end(), f_012 ) - old_faces.begin();
        int index_f_013 = std::find( old_faces.begin(), old_faces.end(), f_013 ) - old_faces.begin();
        
        assert( old_faces[ index_f_012 ] == f_012 );
        assert( old_faces[ index_f_013 ] == f_013 );
        
        int e_0_n = e_0_1;
        int e_n_1 = counter_edges;
        int e_n_2 = counter_edges + 1 + index_f_012;
        int e_n_3 = counter_edges + 1 + index_f_013;
        
        int f_0n2 = f_012;
        int f_n12 = counter_faces + index_f_012;
        int f_0n3 = f_013;
        int f_n13 = counter_faces + index_f_013;
        int f_n23 = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_0n2;
        data_tetrahedron_faces[ t_old ][1] = f_0n3;
        data_tetrahedron_faces[ t_old ][2] = f_023;
        data_tetrahedron_faces[ t_old ][3] = f_n23;
        
        data_tetrahedron_faces[ t_new ][0] = f_n12;
        data_tetrahedron_faces[ t_new ][1] = f_n13;
        data_tetrahedron_faces[ t_new ][2] = f_n23;
        data_tetrahedron_faces[ t_new ][3] = f_123;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_n;
        data_tetrahedron_edges[ t_old ][1] = e_0_2;
        data_tetrahedron_edges[ t_old ][2] = e_0_3;
        data_tetrahedron_edges[ t_old ][3] = e_n_2;
        data_tetrahedron_edges[ t_old ][4] = e_n_3;
        data_tetrahedron_edges[ t_old ][5] = e_2_3;
        
        data_tetrahedron_edges[ t_new ][0] = e_n_1;
        data_tetrahedron_edges[ t_new ][1] = e_n_2;
        data_tetrahedron_edges[ t_new ][2] = e_n_3;
        data_tetrahedron_edges[ t_new ][3] = e_1_2;
        data_tetrahedron_edges[ t_new ][4] = e_1_3;
        data_tetrahedron_edges[ t_new ][5] = e_2_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_n;
        data_tetrahedron_vertices[ t_old ][2] = v_2;
        data_tetrahedron_vertices[ t_old ][3] = v_3;
        
        data_tetrahedron_vertices[ t_new ][0] = v_n;
        data_tetrahedron_vertices[ t_new ][1] = v_1;
        data_tetrahedron_vertices[ t_new ][2] = v_2;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_n_2;
        data_face_edges[ f_new ][1] = e_n_3;
        data_face_edges[ f_new ][2] = e_2_3;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_n;
        data_face_vertices[ f_new ][1] = v_2;
        data_face_vertices[ f_new ][2] = v_3;
        
      } else if( localindex_of_tetrahedron_refinementedge[ ot ] == 1 ) { // 0 2 
        
        assert( v_0 == e_back_vertex && v_2 == e_front_vertex && e == e_0_2 );
        
        /* prepare stuff */ 
        int index_f_012 = std::find( old_faces.begin(), old_faces.end(), f_012 ) - old_faces.begin();
        int index_f_023 = std::find( old_faces.begin(), old_faces.end(), f_023 ) - old_faces.begin();
        
        assert( old_faces[ index_f_012 ] == f_012 );
        assert( old_faces[ index_f_023 ] == f_023 );
        
        int e_0_n = e_0_2;
        int e_n_2 = counter_edges;
        int e_1_n = counter_edges + 1 + index_f_012;
        int e_n_3 = counter_edges + 1 + index_f_023;
        
        int f_01n = f_012;
        int f_1n2 = counter_faces + index_f_012;
        int f_0n3 = f_023;
        int f_2n3 = counter_faces + index_f_023;
        int f_1n3 = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_01n;
        data_tetrahedron_faces[ t_old ][1] = f_013;
        data_tetrahedron_faces[ t_old ][2] = f_0n3;
        data_tetrahedron_faces[ t_old ][3] = f_1n3;
        
        data_tetrahedron_faces[ t_new ][0] = f_1n2;
        data_tetrahedron_faces[ t_new ][1] = f_1n3;
        data_tetrahedron_faces[ t_new ][2] = f_123;
        data_tetrahedron_faces[ t_new ][3] = f_2n3;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_n;
        data_tetrahedron_edges[ t_old ][2] = e_0_3;
        data_tetrahedron_edges[ t_old ][3] = e_1_n;
        data_tetrahedron_edges[ t_old ][4] = e_1_3;
        data_tetrahedron_edges[ t_old ][5] = e_n_3;
        
        data_tetrahedron_edges[ t_new ][0] = e_1_n;
        data_tetrahedron_edges[ t_new ][1] = e_1_2;
        data_tetrahedron_edges[ t_new ][2] = e_1_3;
        data_tetrahedron_edges[ t_new ][3] = e_n_2;
        data_tetrahedron_edges[ t_new ][4] = e_n_3;
        data_tetrahedron_edges[ t_new ][5] = e_2_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_n;
        data_tetrahedron_vertices[ t_old ][3] = v_3;
        
        data_tetrahedron_vertices[ t_new ][0] = v_1;
        data_tetrahedron_vertices[ t_new ][1] = v_n;
        data_tetrahedron_vertices[ t_new ][2] = v_2;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_1_n;
        data_face_edges[ f_new ][1] = e_1_3;
        data_face_edges[ f_new ][2] = e_n_3;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_1;
        data_face_vertices[ f_new ][1] = v_n;
        data_face_vertices[ f_new ][2] = v_3;
        
      } else if( localindex_of_tetrahedron_refinementedge[ ot ] == 2 ) { // 0 3 
        
        assert( v_0 == e_back_vertex && v_3 == e_front_vertex && e == e_0_3 );
        
        /* prepare stuff */ 
        int index_f_013 = std::find( old_faces.begin(), old_faces.end(), f_013 ) - old_faces.begin();
        int index_f_023 = std::find( old_faces.begin(), old_faces.end(), f_023 ) - old_faces.begin();
        
        assert( old_faces[ index_f_013 ] == f_013);
        assert( old_faces[ index_f_023 ] == f_023 );
        
        int e_0_n = e_0_3;
        int e_n_3 = counter_edges;
        int e_1_n = counter_edges + 1 + index_f_013;
        int e_2_n = counter_edges + 1 + index_f_023;
        
        int f_01n = f_013;
        int f_1n3 = counter_faces + index_f_013;
        int f_02n = f_023;
        int f_2n3 = counter_faces + index_f_023;
        int f_12n = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_012;
        data_tetrahedron_faces[ t_old ][1] = f_01n;
        data_tetrahedron_faces[ t_old ][2] = f_02n;
        data_tetrahedron_faces[ t_old ][3] = f_12n;
        
        data_tetrahedron_faces[ t_new ][0] = f_12n;
        data_tetrahedron_faces[ t_new ][1] = f_123;
        data_tetrahedron_faces[ t_new ][2] = f_1n3;
        data_tetrahedron_faces[ t_new ][3] = f_2n3;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_2;
        data_tetrahedron_edges[ t_old ][2] = e_0_n;
        data_tetrahedron_edges[ t_old ][3] = e_1_2;
        data_tetrahedron_edges[ t_old ][4] = e_1_n;
        data_tetrahedron_edges[ t_old ][5] = e_2_n;
        
        data_tetrahedron_edges[ t_new ][0] = e_1_2;
        data_tetrahedron_edges[ t_new ][1] = e_1_n;
        data_tetrahedron_edges[ t_new ][2] = e_1_3;
        data_tetrahedron_edges[ t_new ][3] = e_2_n;
        data_tetrahedron_edges[ t_new ][4] = e_2_3;
        data_tetrahedron_edges[ t_new ][5] = e_n_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_2;
        data_tetrahedron_vertices[ t_old ][3] = v_n;
        
        data_tetrahedron_vertices[ t_new ][0] = v_1;
        data_tetrahedron_vertices[ t_new ][1] = v_2;
        data_tetrahedron_vertices[ t_new ][2] = v_n;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_1_2;
        data_face_edges[ f_new ][1] = e_1_n;
        data_face_edges[ f_new ][2] = e_2_n;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_1;
        data_face_vertices[ f_new ][1] = v_2;
        data_face_vertices[ f_new ][2] = v_n;
    
      } else if( localindex_of_tetrahedron_refinementedge[ ot ] == 3 ) { // 1 2 
        
        assert( v_1 == e_back_vertex && v_2 == e_front_vertex && e == e_1_2 );
        
        /* prepare stuff */ 
        int index_f_012 = std::find( old_faces.begin(), old_faces.end(), f_012 ) - old_faces.begin();
        int index_f_123 = std::find( old_faces.begin(), old_faces.end(), f_123 ) - old_faces.begin();
        
        assert( old_faces[ index_f_012 ] == f_012 );
        assert( old_faces[ index_f_123 ] == f_123 );
        
        int e_1_n = e_1_2;
        int e_n_2 = counter_edges;
        int e_0_n = counter_edges + 1 + index_f_012;
        int e_n_3 = counter_edges + 1 + index_f_123;
        
        int f_01n = f_012;
        int f_0n2 = counter_faces + index_f_012;
        int f_1n3 = f_123;
        int f_n23 = counter_faces + index_f_123;
        int f_0n3 = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_01n;
        data_tetrahedron_faces[ t_old ][1] = f_013;
        data_tetrahedron_faces[ t_old ][2] = f_0n3;
        data_tetrahedron_faces[ t_old ][3] = f_1n3;
        
        data_tetrahedron_faces[ t_new ][0] = f_0n2;
        data_tetrahedron_faces[ t_new ][1] = f_0n3;
        data_tetrahedron_faces[ t_new ][2] = f_023;
        data_tetrahedron_faces[ t_new ][3] = f_n23;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_n;
        data_tetrahedron_edges[ t_old ][2] = e_0_3;
        data_tetrahedron_edges[ t_old ][3] = e_1_n;
        data_tetrahedron_edges[ t_old ][4] = e_1_3;
        data_tetrahedron_edges[ t_old ][5] = e_n_3;
        
        data_tetrahedron_edges[ t_new ][0] = e_0_n;
        data_tetrahedron_edges[ t_new ][1] = e_0_2;
        data_tetrahedron_edges[ t_new ][2] = e_0_3;
        data_tetrahedron_edges[ t_new ][3] = e_n_2;
        data_tetrahedron_edges[ t_new ][4] = e_n_3;
        data_tetrahedron_edges[ t_new ][5] = e_2_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_n;
        data_tetrahedron_vertices[ t_old ][3] = v_3;
        
        data_tetrahedron_vertices[ t_new ][0] = v_0;
        data_tetrahedron_vertices[ t_new ][1] = v_n;
        data_tetrahedron_vertices[ t_new ][2] = v_2;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_0_n;
        data_face_edges[ f_new ][1] = e_0_3;
        data_face_edges[ f_new ][2] = e_n_3;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_0;
        data_face_vertices[ f_new ][1] = v_n;
        data_face_vertices[ f_new ][2] = v_3;
    
      } else if( localindex_of_tetrahedron_refinementedge[ ot ] == 4 ) { // 1 3 
        
        assert( v_1 == e_back_vertex && v_3 == e_front_vertex && e == e_1_3 );
        
        /* prepare stuff */ 
        int index_f_013 = std::find( old_faces.begin(), old_faces.end(), f_013 ) - old_faces.begin();
        int index_f_123 = std::find( old_faces.begin(), old_faces.end(), f_123 ) - old_faces.begin();
        
        assert( old_faces[ index_f_013 ] == f_013 );
        assert( old_faces[ index_f_123 ] == f_123 );
        
        int e_1_n = e_1_3;
        int e_n_3 = counter_edges;
        int e_0_n = counter_edges + 1 + index_f_013;
        int e_2_n = counter_edges + 1 + index_f_123;
        
        int f_01n = f_013;
        int f_0n3 = counter_faces + index_f_013;
        int f_12n = f_123;
        int f_2n3 = counter_faces + index_f_123;
        int f_02n = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_012;
        data_tetrahedron_faces[ t_old ][1] = f_01n;
        data_tetrahedron_faces[ t_old ][2] = f_02n;
        data_tetrahedron_faces[ t_old ][3] = f_12n;
        
        data_tetrahedron_faces[ t_new ][0] = f_02n;
        data_tetrahedron_faces[ t_new ][1] = f_023;
        data_tetrahedron_faces[ t_new ][2] = f_0n3;
        data_tetrahedron_faces[ t_new ][3] = f_2n3;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_2;
        data_tetrahedron_edges[ t_old ][2] = e_0_n;
        data_tetrahedron_edges[ t_old ][3] = e_1_2;
        data_tetrahedron_edges[ t_old ][4] = e_1_n;
        data_tetrahedron_edges[ t_old ][5] = e_2_n;
        
        data_tetrahedron_edges[ t_new ][0] = e_0_2;
        data_tetrahedron_edges[ t_new ][1] = e_0_n;
        data_tetrahedron_edges[ t_new ][2] = e_0_3;
        data_tetrahedron_edges[ t_new ][3] = e_2_n;
        data_tetrahedron_edges[ t_new ][4] = e_2_3;
        data_tetrahedron_edges[ t_new ][5] = e_n_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_2;
        data_tetrahedron_vertices[ t_old ][3] = v_n;
        
        data_tetrahedron_vertices[ t_new ][0] = v_0;
        data_tetrahedron_vertices[ t_new ][1] = v_2;
        data_tetrahedron_vertices[ t_new ][2] = v_n;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_0_2;
        data_face_edges[ f_new ][1] = e_0_n;
        data_face_edges[ f_new ][2] = e_2_n;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_0;
        data_face_vertices[ f_new ][1] = v_2;
        data_face_vertices[ f_new ][2] = v_n;
    
      } else if( localindex_of_tetrahedron_refinementedge[ ot ] == 5 ) { // 2 3 
        
        assert( v_2 == e_back_vertex && v_3 == e_front_vertex && e == e_2_3 );
        
        /* prepare stuff */ 
        int index_f_023 = std::find( old_faces.begin(), old_faces.end(), f_023 ) - old_faces.begin();
        int index_f_123 = std::find( old_faces.begin(), old_faces.end(), f_123 ) - old_faces.begin();
        
        assert( old_faces[ index_f_023 ] == f_023 );
        assert( old_faces[ index_f_123 ] == f_123 );
        
        int e_2_n = e_2_3;
        int e_n_3 = counter_edges;
        int e_0_n = counter_edges + 1 + index_f_023;
        int e_1_n = counter_edges + 1 + index_f_123;
        
        int f_02n = f_023;
        int f_0n3 = counter_faces + index_f_023;
        int f_12n = f_123;
        int f_1n3 = counter_faces + index_f_123;
        int f_01n = counter_faces + old_faces.size() + ot;
        
        /* Tetrahedron -> Faces */
        data_tetrahedron_faces[ t_old ][0] = f_012;
        data_tetrahedron_faces[ t_old ][1] = f_01n;
        data_tetrahedron_faces[ t_old ][2] = f_02n;
        data_tetrahedron_faces[ t_old ][3] = f_12n;
        
        data_tetrahedron_faces[ t_new ][0] = f_01n;
        data_tetrahedron_faces[ t_new ][1] = f_013;
        data_tetrahedron_faces[ t_new ][2] = f_0n3;
        data_tetrahedron_faces[ t_new ][3] = f_1n3;
        
        /* Tetrahedron -> Edges */
        data_tetrahedron_edges[ t_old ][0] = e_0_1;
        data_tetrahedron_edges[ t_old ][1] = e_0_2;
        data_tetrahedron_edges[ t_old ][2] = e_0_n;
        data_tetrahedron_edges[ t_old ][3] = e_1_2;
        data_tetrahedron_edges[ t_old ][4] = e_1_n;
        data_tetrahedron_edges[ t_old ][5] = e_2_n;
        
        data_tetrahedron_edges[ t_new ][0] = e_0_1;
        data_tetrahedron_edges[ t_new ][1] = e_0_n;
        data_tetrahedron_edges[ t_new ][2] = e_0_3;
        data_tetrahedron_edges[ t_new ][3] = e_1_n;
        data_tetrahedron_edges[ t_new ][4] = e_1_3;
        data_tetrahedron_edges[ t_new ][5] = e_n_3;
        
        /* Tetrahedron -> Vertices */
        data_tetrahedron_vertices[ t_old ][0] = v_0;
        data_tetrahedron_vertices[ t_old ][1] = v_1;
        data_tetrahedron_vertices[ t_old ][2] = v_2;
        data_tetrahedron_vertices[ t_old ][3] = v_n;
        
        data_tetrahedron_vertices[ t_new ][0] = v_0;
        data_tetrahedron_vertices[ t_new ][1] = v_1;
        data_tetrahedron_vertices[ t_new ][2] = v_n;
        data_tetrahedron_vertices[ t_new ][3] = v_3;
        
        /* Face -> Edges */
        data_face_edges[ f_new ][0] = e_0_1;
        data_face_edges[ f_new ][1] = e_0_n;
        data_face_edges[ f_new ][2] = e_1_n;
        
        /* Face -> Vertices */
        data_face_vertices[ f_new ][0] = v_0;
        data_face_vertices[ f_new ][1] = v_1;
        data_face_vertices[ f_new ][2] = v_n;
    
      } else {
        
        unreachable();
        
      }
       
      
    }
      
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    /***********************************************/
    /*****                                    ******/
    /*****   Delete the adjaceny links        ******/
    /*****   in the macropatch                ******/
    /*****   regarding micropatch simplices   ******/
    /*****                                    ******/
    /***********************************************/
    
    for( int v : micropatch_vertices ) data_vertex_firstparent_edge[v] = nullindex;
    for( int v : micropatch_vertices ) data_vertex_firstparent_face[v] = nullindex;
    for( int v : micropatch_vertices ) data_vertex_firstparent_tetrahedron[v] = nullindex;
    for( int e : micropatch_edges ) data_edge_firstparent_face[e] = nullindex;
    for( int e : micropatch_edges ) data_edge_firstparent_tetrahedron[e] = nullindex;
    for( int f : micropatch_faces ) data_face_firstparent_tetrahedron[f] = nullindex;
    
    for( int e : macropatch_edges ) 
    for( int vi = 0; vi < 2; vi++ ) 
        if( micropatch_vertices.find( data_edge_vertices[e][vi] ) != micropatch_vertices.end() )
            data_edge_nextparents_of_vertices[e][vi] = nullindex;
        
    for( int f : macropatch_faces ) 
    for( int vi = 0; vi < 3; vi++ ) 
        if( micropatch_vertices.find( data_face_vertices[f][vi] ) != micropatch_vertices.end() )
            data_face_nextparents_of_vertices[f][vi] = nullindex;
        
    for( int f : macropatch_faces ) 
    for( int ei = 0; ei < 3; ei++ ) 
        if( micropatch_edges.find( data_face_edges[f][ei] ) != micropatch_edges.end() )
            data_face_nextparents_of_edges[f][ei] = nullindex;
        
         
    for( int t : macropatch_tetrahedra ) 
    for( int vi = 0; vi < 4; vi++ ) 
        if( micropatch_vertices.find( data_tetrahedron_vertices[t][vi] ) != micropatch_vertices.end() )
            data_tetrahedron_nextparents_of_vertices[t][vi] = nullindex;
        
    for( int t : macropatch_tetrahedra ) 
    for( int ei = 0; ei < 6; ei++ ) 
        if( micropatch_edges.find( data_tetrahedron_edges[t][ei] ) != micropatch_edges.end() )
            data_tetrahedron_nextparents_of_edges[t][ei] = nullindex;
       
    for( int t : macropatch_tetrahedra ) 
    for( int fi = 0; fi < 4; fi++ ) 
        if( micropatch_faces.find( data_tetrahedron_faces[t][fi] ) != micropatch_faces.end() )
            data_tetrahedron_nextparents_of_faces[t][fi] = nullindex;
    
    
        
        
        
        
    
    
    /***********************************************/
    /*****                                    ******/
    /*****   Rebuild the adjaceny links       ******/
    /*****   in the macropatch                ******/
    /*****   regarding micropatch simplices   ******/
    /*****                                    ******/
    /***********************************************/
    
    for( int e : macropatch_edges )
    for( int vi = 0; vi < 2; vi++ )
    {
        
        int v = data_edge_vertices[e][vi];
        
        if( micropatch_vertices.find( v ) != micropatch_vertices.end() ){
            
            int fp = data_vertex_firstparent_edge[v]; 
            
            data_vertex_firstparent_edge[v] = e;
            
            data_edge_nextparents_of_vertices[e][vi] = fp;
            
        }
        
    }
    
    
    for( int f : macropatch_faces )
    for( int vi = 0; vi < 3; vi++ )
    {
        
        int v = data_face_vertices[f][vi];
        
        if( micropatch_vertices.find( v ) != micropatch_vertices.end() ){
            
            int fp = data_vertex_firstparent_face[v]; 
            
            data_vertex_firstparent_face[v] = f;
            
            data_face_nextparents_of_vertices[f][vi] = fp;
            
        }
        
    }
    
    for( int f : macropatch_faces )
    for( int ei = 0; ei < 3; ei++ )
    {
        
        int e = data_face_edges[f][ei];
        
        if( micropatch_edges.find( e ) != micropatch_edges.end() ){
            
            int fp = data_edge_firstparent_face[e]; 
            
            data_edge_firstparent_face[e] = f;
            
            data_face_nextparents_of_edges[f][ei] = fp;
            
        }
        
    }
    
    
    
    
    for( int t : macropatch_tetrahedra )
    for( int vi = 0; vi < 4; vi++ )
    {
        
        int v = data_tetrahedron_vertices[t][vi];
        
        if( micropatch_vertices.find( v ) != micropatch_vertices.end() ){
            
            int fp = data_vertex_firstparent_tetrahedron[v]; 
            
            data_vertex_firstparent_tetrahedron[v] = t;
            
            data_tetrahedron_nextparents_of_vertices[t][vi] = fp;
            
        }
        
    }
    
    for( int t : macropatch_tetrahedra )
    for( int ei = 0; ei < 6; ei++ )
    {
        
        int e = data_tetrahedron_edges[t][ei];
        
        if( micropatch_edges.find( e ) != micropatch_edges.end() ){
            
            int fp = data_edge_firstparent_tetrahedron[e]; 
            
            data_edge_firstparent_tetrahedron[e] = t;
            
            data_tetrahedron_nextparents_of_edges[t][ei] = fp;
            
        }
        
    }
    
    for( int t : macropatch_tetrahedra )
    for( int fi = 0; fi < 4; fi++ )
    {
        
        int f = data_tetrahedron_faces[t][fi];
        
        if( micropatch_faces.find( f ) != micropatch_faces.end() ){
            
            int fp = data_face_firstparent_tetrahedron[f]; 
            
            data_face_firstparent_tetrahedron[f] = t;
            
            data_tetrahedron_nextparents_of_faces[t][fi] = fp;
            
        }
        
    }
    
    
    /**************************/
    /*     ADJACENCY LINKS    */
    /**************************/
    
    // Run over the simplices in the macropatch 
    // and reassemble the links corresponding 
    // to simplices in the micropatch
    // or newly created. 
    
    
    /**************************************************/
    
    
    /*
     *  UPDATE COUNTERS 
     */
    
    counter_tetrahedra += old_tetrahedra.size();
    counter_faces      += old_faces.size() + old_tetrahedra.size();
    counter_edges      += 1 + old_faces.size();
    counter_vertices   += 1;
    
    
    
    
    
    /* Done */
    
    // check();
    
}

/*
 * * * * * VERTEX BISECTION
 */







/*
void MeshSimplicial3D::longest_edge_bisection( std::vector<int> edges )
{
    check();
    
    const int old_vertex_count = counter_vertices;
    
    
    // 0. check the input
    
    for( int& e : edges )
        assert( 0 <= e && e < counter_edges );
    
    
    // 1. create stack for the edges to be bisected, and fill in first batch 
    
    std::stack<int> todostack;
    
    for( int& e : edges )
        todostack.push( e );
        
    // 2. conduct the main loop of the refinement algorithm 
    
    while( ! todostack.empty() )
    {
        
        std::cout << todostack.size() << " ";
        
        // as long the stack is not empty,
        // pick the top edge and make the following case distinction
        // a) the edge index belongs to an edge already bisected and can be ignored 
        // b) if the top edge is bisection edge to all its neighbors, bisect and pop
        // c) else, push the necessary edge of each parent simplex
        
        int e = todostack.top();
        
        // to check whether e belongs to an edge that has already been bisected,
        // we check whether one of the vertices belongs to the new vertices 
        
        if( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count ) {
            
            todostack.pop();
        
        } else {
            
            // run over neighbor tetrahedra and check for necessary edges 
            
            for( int t = get_edge_firstparent_tetrahedron( e ); t != nullindex; t = get_edge_nextparent_tetrahedron( e, t ) )
            {
                
                int e0 = get_tetrahedron_edge( t, 0 ); Float l0 = get_edge_length( e0 );
                int e1 = get_tetrahedron_edge( t, 1 ); Float l1 = get_edge_length( e1 );
                int e2 = get_tetrahedron_edge( t, 2 ); Float l2 = get_edge_length( e2 );
                int e3 = get_tetrahedron_edge( t, 3 ); Float l3 = get_edge_length( e3 );
                int e4 = get_tetrahedron_edge( t, 4 ); Float l4 = get_edge_length( e4 );
                int e5 = get_tetrahedron_edge( t, 5 ); Float l5 = get_edge_length( e5 );
                
                assert( e == e0 || e == e1 || e == e2 || e == e3 || e == e4 || e == e5 );
                
                int ne = nullindex;
                if( l0 >= l1 && l0 >= l2 && l0 >= l3 && l0 >= l4 && l0 >= l5 ) ne = e0;
                if( l1 >= l0 && l1 >= l2 && l1 >= l3 && l1 >= l4 && l1 >= l5 ) ne = e1;
                if( l2 >= l0 && l2 >= l1 && l2 >= l3 && l2 >= l4 && l2 >= l5 ) ne = e2;
                if( l3 >= l0 && l3 >= l1 && l3 >= l2 && l3 >= l4 && l3 >= l5 ) ne = e3;
                if( l4 >= l0 && l4 >= l1 && l4 >= l2 && l4 >= l3 && l4 >= l5 ) ne = e4;
                if( l5 >= l0 && l5 >= l1 && l5 >= l2 && l5 >= l3 && l5 >= l4 ) ne = e5;
                assert( ne != nullindex );
                
                if( ne != e )
                    todostack.push( ne );
                
            }
            
            // if top edge is still the same, bisect
            
            if( e == todostack.top() )
            {
                todostack.pop();
                bisect_edge( e );
            }
            
        }
        
    }
    
    // fiinished!
    
    check();
}*/





// void MeshSimplicial3D::longest_edge_bisection( std::vector<int> edges )
// {
//     check();
    
//     const int old_vertex_count = counter_vertices;
    
    
//     // 0. check the input
    
//     for( int& e : edges )
//         assert( 0 <= e && e < counter_edges );
    
    
//     // 1. create stack for the edges to be bisected, and fill in first batch 
    
//     std::list<int> todostack;
    
//     for( int& e : edges )
//         todostack.push_back( e ); // put e on top either by inserting or pulling it up!
        
    
//     // 2. conduct the main loop of the refinement algorithm 
    
//     while( ! todostack.empty() )
//     {
        
//         //LOG << todostack.back() << space << todostack.size();
        
//         // as long the stack is not empty,
//         // pick the top edge and make the following case distinction
//         // a) if the top edge is longer than its neighbors, bisect and pop
//         // b) else, push the longest edge of each parent simplex
        
//         int e = todostack.back();
    
//         // to check whether e belongs to an edge that has already been bisected,
//         // we check whether one of the vertices belongs to the new vertices 
        
//         // if( get_edge_vertex( e, 0 ) >= old_vertex_count || get_edge_vertex(e,1) >= old_vertex_count ) {
            
//         //     todostack.pop();
        
//         // } else 
//         {
            
//             Float length_e = get_edge_length( e );
                
//             // run over neighbor tetrahedra and check for longer edges 
            
//             bool compatibly_divisible = true;
            
//             for( int t = get_edge_firstparent_tetrahedron( e ); t != nullindex; t = get_edge_nextparent_tetrahedron( e, t ) )
//             for( int ei = 0; ei < 6; ei++ )
//             {
//                 int other_e = get_tetrahedron_edge( t, ei );
                
//                 if( e != other_e && get_edge_length( other_e ) > length_e ) {
//                     compatibly_divisible = false;
//                     todostack.push_back( other_e );
//                 }
//             }
            
//             // if top edge is still the same, bisect 
        
//             if( compatibly_divisible )
//             {
//                 assert( e == todostack.back() );
//                 todostack.remove( e );
//                 bisect_edge( e );
//                 assert( std::find( todostack.begin(), todostack.end(), e ) == todostack.end() );
//             } else 
//                 assert( e != todostack.back() );

//         }
           
//     }
    
    
//     // fiinished!
    
//     check();
// }







/*
 * * * * * LONGEST EDGE VERTEX BISECTION
 */

void MeshSimplicial3D::longest_edge_bisection_recursive( const std::vector<int>& edges )
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


void MeshSimplicial3D::longest_edge_bisection_recursive( const int e )
{
    assert( 0 <= e && e < counter_edges );

    int longest_edge = nullindex;
        
    do{

        longest_edge = e;
        
        for( int t = get_edge_firstparent_tetrahedron( e ); t != nullindex; t = get_edge_nextparent_tetrahedron( e, t ) )
        for( int ei = 0; ei < 6; ei++ )
        {
            int other_e = get_tetrahedron_edge( t, ei );
            
            if( e != other_e && get_edge_length( other_e ) > get_edge_length( e ) )
                longest_edge = other_e;
        }

        if( longest_edge != e )
            longest_edge_bisection_recursive( longest_edge );

    } while ( longest_edge != e );

    bisect_edge( e );

}


