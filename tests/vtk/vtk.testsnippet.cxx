
#include "../../utility/utility.hpp"
#include "../../operators/floatvector.hpp"
#include "../../mesh/coordinates.hpp"


    
    
inline void internal_print( const Mesh& M, std::string meshname, std::string filename = "" )
{
    
    std::string basename = ( filename == "" ) ? getbasename(__FILE__) : filename;
    
    fstream fs( experimentfile( basename ), std::fstream::out );
    
    VTKWriter vtk( M, fs, meshname );
    // vtk.writeCoordinateBlock();
    // vtk.writeTopDimensionalCells();
    
    {
        FloatVector V( M.count_simplices(0), 
                        [&M](int i)->Float{
                        FloatVector point = M.getcoordinates().getvectorclone(i);
                        // return i / (Float) M.count_simplices(0);
                        Float x = point[0];
                        Float y = ( ( point.getdimension() >= 2 ) ? point[1] : 1. );
                        Float z = ( ( point.getdimension() >= 3 ) ? point[2] : 1. );
                        return std::sin( Constants::pi * x ) * y + z;
                    });
        
        vtk.writeVertexScalarData( V, "testing_scalar_data", 1.0 );
    }

    int inner_dim = M.getinnerdimension();
    
    vtk.writeCellScalarData( FloatVector(M.count_simplices(inner_dim), 2.5 ), "cell_scalar_data" );
    
    vtk.writeCellVectorData( 
        FloatVector( M.count_simplices(inner_dim),  1.5 ),
        FloatVector( M.count_simplices(inner_dim),  2.5 ),
        FloatVector( M.count_simplices(inner_dim), -3.5 ),
        "cell_vector_data"
    );

    fs.close();
}
