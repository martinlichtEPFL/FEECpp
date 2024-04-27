#ifndef INCLUDEGUARD_MESH_IO_SIMPLICIAL3D_HPP
#define INCLUDEGUARD_MESH_IO_SIMPLICIAL3D_HPP


#include <istream>
#include <ostream>


#include "../basic.hpp"
#include "mesh.simplicial3D.hpp"


/*******************
****  
****  Input / Output of coordinate objects 
****  
*******************/


void writeMeshSimplicial3D( std::ostream& out, const MeshSimplicial3D& mesh, bool sugar = false );

MeshSimplicial3D readMeshSimplicial3D( std::istream& in );

void writeMeshSimplicial3D( const char* filename, const MeshSimplicial3D& mesh, bool sugar = false );

MeshSimplicial3D readMeshSimplicial3D( const char* filename );



#endif
