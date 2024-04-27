
#include "vtkwriter.hpp"
#include "../combinatorics/generateindexmaps.hpp"


VTKWriter::VTKWriter( const Mesh& m, std::ostream& os, const std::string& name )
: VTKWriter( m, os, name, [&](int i) -> Float { return 0.; } )
{}

VTKWriter::VTKWriter( const Mesh& m, std::ostream& os, const std::string& name, const FloatVector& z )
: VTKWriter( m, os, name, [&](int i) -> Float { return z[i]; } )
{
    assert( z.isfinite() );
}

VTKWriter::VTKWriter( const Mesh& m, std::ostream& os, const std::string& name, const std::function<Float(int)>& func_z )
: mesh(m), os(os), current_stage(Stage::nothing)
{
    mesh.check();
    
    assert( mesh.getouterdimension() == 1 || mesh.getouterdimension() == 2 || mesh.getouterdimension() == 3 );
    assert( mesh.getinnerdimension() == 1 || mesh.getinnerdimension() == 2 || mesh.getinnerdimension() == 3 );
    
    assert( mesh.has_dimension_counted( mesh.getinnerdimension() ) );
    assert( mesh.has_dimension_counted(0) );
    assert( mesh.has_subsimplices_listed( mesh.getinnerdimension(), 0 ) );
    
    writePreamble( name );

    assert( current_stage == Stage::preamble );

    writeCoordinateBlock( func_z );

    assert( current_stage == Stage::coordinate );

    writeTopDimensionalCells();
    
    assert( current_stage == Stage::cells );

}


VTKWriter VTKWriter::writePreamble( const std::string& name )
{
    assert( current_stage == Stage::nothing );
    current_stage = Stage::preamble;
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    
    os << "# vtk DataFile Version 3.0" << nl;
    os << name << nl;
    os << "ASCII" << nl;
    os << "DATASET UNSTRUCTURED_GRID" << nl;
    
    return *this;
}



VTKWriter VTKWriter::writeCoordinateBlock()
{
    writeCoordinateBlock( [](int)->Float { return 0.; } );
    return *this;
}
        
VTKWriter VTKWriter::writeCoordinateBlock( const FloatVector& z )
{
    assert( mesh.count_simplices(0) == z.getdimension() );
    assert( z.isfinite() );
    writeCoordinateBlock( [&](int i)->Float { return z.at(i); } );
    return *this;
}
        
        
VTKWriter VTKWriter::writeCoordinateBlock( const std::function<Float(int)>& func_z )
{
    assert( current_stage == Stage::preamble );
    current_stage = Stage::coordinate;
    
    os << "POINTS " << mesh.count_simplices(0) << " double" << nl;
    
    for( int v = 0; v < mesh.count_simplices(0); v++ )
        if( mesh.getouterdimension() == 1 ) {
            os << mesh.getcoordinates().getdata(v,0)
               << space
               << 0.0 
               << space 
               << func_z(v)
               << nl;
        } else if( mesh.getouterdimension() == 2 ) {
            os << mesh.getcoordinates().getdata(v,0)
               << space
               << mesh.getcoordinates().getdata(v,1) 
               << space 
               << func_z(v)
               << nl;
        } else if( mesh.getouterdimension() == 3 ) {
            os << mesh.getcoordinates().getdata(v,0)
               << space
               << mesh.getcoordinates().getdata(v,1) 
               << space 
               << mesh.getcoordinates().getdata(v,2) 
               << nl;
        } else {
            unreachable();
        }
        
    
    os << nl;

    return *this;
}
        
        
VTKWriter VTKWriter::writeTopDimensionalCells()
{
    assert( current_stage == Stage::coordinate);
    current_stage = Stage::cells;
    
    int topdim = mesh.getinnerdimension();
    
    assert( 1 <= topdim and topdim <= 3 );
    
    os << "CELLS " 
       << mesh.count_simplices(topdim)
       << space 
       << (topdim+1) * mesh.count_simplices(topdim) + mesh.count_simplices(topdim)
       << nl;
    
    for( int S = 0; S < mesh.count_simplices(topdim); S++ ) {
        os << topdim+1; 
        const auto list_of_vertices = mesh.getsubsimplices(topdim,0,S);
        for( int i = 0; i <= topdim; i++ ) os << space << list_of_vertices[i];
        os << nl;
    }
    
    os << nl;
    
    os << "CELL_TYPES" << space << mesh.count_simplices(topdim) << nl;
    
    if( topdim == 1 ) for( int S = 0; S < mesh.count_simplices(1); S++ ) os <<  3 << nl;
    if( topdim == 2 ) for( int S = 0; S < mesh.count_simplices(2); S++ ) os <<  5 << nl;
    if( topdim == 3 ) for( int S = 0; S < mesh.count_simplices(3); S++ ) os << 10 << nl;
    
    os << nl;

    return *this;
}

















VTKWriter VTKWriter::writeVertexScalarData( const std::function<Float(int)>& datafunction, const std::string& name, Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::fielddata );
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    if( current_stage != Stage::fielddata ){
        os << "POINT_DATA " << mesh.count_simplices(0) << nl << nl;
        current_stage = Stage::fielddata;
    }
    
    os << "SCALARS " << name << " double 1" << nl;
    // SCALARS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1)
    os << "LOOKUP_TABLE default" << nl;
    
    for( int v = 0; v < mesh.count_simplices(0); v++ )
        os << scaling * datafunction(v) << nl;
    
    os << nl;

    return *this;
}

VTKWriter VTKWriter::writeCellScalarData( const std::function<Float(int)>& datafunction, const std::string& name, Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::fielddata );
    
    const int topdim = mesh.getinnerdimension();
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    if( current_stage != Stage::fielddata ){
        os << "CELL_DATA " << mesh.count_simplices(topdim) << nl << nl;
        current_stage = Stage::fielddata;
    }
    
    os << "SCALARS " << name << " double 1" << nl;
    // VECTORS (name_of_data) Datentyp(=float) AnzahlKomponenten(=1) 
    os << "LOOKUP_TABLE default" << nl;
    
    for( int c = 0; c < mesh.count_simplices(topdim); c++ )
        os << scaling * datafunction(c) << nl;
    
    os << nl;

    return *this;
}
















VTKWriter VTKWriter::writeVertexScalarData( const std::function<Float(const FloatVector&)>& function, const std::string& name, Float scaling )
{
    auto datafunction = [&](int v) -> Float { 
        assert( 0 <= v and v < this->mesh.count_simplices(0) );
        return function( mesh.get_midpoint(0,v) ); 
    };

    writeVertexScalarData( datafunction, name, scaling );

    return *this;
}

VTKWriter VTKWriter::writeVertexScalarData( const std::function<FloatVector(const FloatVector&)>& function, const std::string& name, Float scaling )
{
    auto datafunction = [&](int v) -> Float { 
        assert( 0 <= v and v < this->mesh.count_simplices(0) );
        FloatVector evaluation = function( mesh.get_midpoint(0,v) ); 
        assert( evaluation.getdimension() == 1 );
        return evaluation[0];
    };

    writeVertexScalarData( datafunction, name, scaling );

    return *this;
}

VTKWriter VTKWriter::writeVertexScalarData( const FloatVector& data, const std::string& name, Float scaling )
{
    assert( data.getdimension() == mesh.count_simplices(0) );
    assert( data.isfinite() );

    auto datafunction = [&](int v) -> Float { 
        assert( 0 <= v and v < this->mesh.count_simplices(0) );
        return data.at(v); 
    };

    writeVertexScalarData( datafunction, name, scaling );
    
    return *this;
}






VTKWriter VTKWriter::writeCellScalarData( const std::function<Float(const FloatVector&)>& function, const std::string& name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();

    auto datafunction = [&](int c) -> Float { 
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        return function( mesh.get_midpoint(topdim,c) ); 
    };

    writeCellScalarData( datafunction, name, scaling );

    return *this;
}

VTKWriter VTKWriter::writeCellScalarData( const std::function<FloatVector(const FloatVector&)>& function, const std::string& name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();

    auto datafunction = [&](int c) -> Float { 
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        FloatVector evaluation = function( mesh.get_midpoint(topdim,c) ); 
        assert( evaluation.getdimension() == 1 );
        return evaluation[0];
    };

    writeCellScalarData( datafunction, name, scaling );

    return *this;
}

VTKWriter VTKWriter::writeCellScalarData( const FloatVector& data, const std::string& name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();

    assert( data.getdimension() == mesh.count_simplices(topdim) );
    assert( data.isfinite() );

    auto datafunction = [&](int c) -> Float { 
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        return data.at(c); 
    };

    writeCellScalarData( datafunction, name, scaling );
    
    return *this;
}

VTKWriter VTKWriter::writeCellScalarData_barycentricvolumes( const FloatVector& v, const std::string& name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();
    
    assert( v.getdimension() == mesh.count_simplices(topdim) * (mesh.getinnerdimension()+1) );
    assert( v.isfinite() );

    auto sigmas = generateSigmas( IndexRange(1,topdim), IndexRange(0,topdim) );

    std::vector<Float> signs( sigmas.size(), notanumber );

    // If zero is in the range of a sigma and index i is missing,
    // then we first expand 0 to the sum of its complements (one sign).
    // Only the summand with index i remains,
    // and we order into the remaining indices (i-1 swaps).
    // In total, if index i missing, we get i signs
    // 
    // That is the same sign as sorting sigma^c and sigma
    // 
    for( int i = 0; i <= topdim; i++ )
        signs.at(i) = sign_of_rho_sigma( sigmas[i] ); 

    auto datafunction = [&](int c) -> Float {
    
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        
        Float coefficient = 0.;
        
        for( int i = 0; i <= mesh.getinnerdimension(); i++ ) 
            coefficient += signs.at(i) * v.at( c * (topdim+1) + i ); 
        
        // The barycentric form d\lambda_1 \wedge \dots \wedge d\lambda_n
        // equals the euclidean form dx_1 \wedge \dots \wedge dx_n,
        // up to multiplication by the determinant of the transformation jacobian.
        // That Jacobian is measure(T)/(n+1)!
        // 
        Float factor_bary_to_euclid = Determinant( mesh.getTransformationJacobian(topdim,c) );
        //mesh.getMeasure(topdim,c) / factorial_numerical(topdim);
        Float direction = scaling * factor_bary_to_euclid * coefficient; 

        return direction;
    };
    
    writeCellScalarData( datafunction, name, scaling );
    
    return *this;
}













VTKWriter VTKWriter::writeCellVectorData( const std::function<FloatVector(int)>& datafunction, const std::string& name, Float scaling )
{
    assert( current_stage >= Stage::cells );
    assert( current_stage <= Stage::fielddata );
    
    assert( 0 < name.size() and name.size() <= 256 );
    assert( name.find('\n') == std::string::npos );
    assert( count_white_space( name ) == 0 );
    
    const int topdim = mesh.getinnerdimension();
    
    if( current_stage != Stage::fielddata ){
        os << "CELL_DATA " << mesh.count_simplices(topdim) << nl << nl;
        current_stage = Stage::fielddata;
    }
    
    os << "VECTORS " << name << " double" << nl;
    // VECTORS (name_of_data) Datentyp(=float) 
    
    for( int c = 0; c < mesh.count_simplices(topdim); c++ ) {

        auto data_item = scaling * datafunction(c);

        os << data_item.at(0) << space 
           << ( mesh.getouterdimension() >= 2 ? data_item.at(1) : 0. ) << space 
           << ( mesh.getouterdimension() >= 3 ? data_item.at(2) : 0. ) << nl;
    }
    
    os << nl;

    return *this;
}


VTKWriter VTKWriter::writeCellVectorData( 
    const FloatVector& datax,
    const FloatVector& datay,
    const FloatVector& dataz,
    const std::string& name, 
    Float scaling )
{
    const int topdim = mesh.getinnerdimension();
    
    assert( datax.getdimension() == datay.getdimension() );
    assert( datax.getdimension() == dataz.getdimension() );
    assert( datax.getdimension() == mesh.count_simplices(topdim) );

    assert( datax.isfinite() );
    assert( datay.isfinite() );
    assert( dataz.isfinite() );

    
    auto datafunction = [&](int c) -> FloatVector { 
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        return FloatVector({datax[c],datay[c],dataz[c]}); 
    };

    writeCellVectorData( datafunction, name, scaling );
    
    return *this;
}



VTKWriter VTKWriter::writeCellVectorData( const std::function<FloatVector(const FloatVector&)>& function, const std::string& name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();

    auto datafunction = [&](int c) -> FloatVector { 
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        return function( mesh.get_midpoint(topdim,c) ); 
    };
    
    writeCellVectorData( datafunction, name, scaling );

    return *this;
}















VTKWriter VTKWriter::writeCellVectorData_barycentricgradients( const FloatVector& v, const std::string& name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();
    
    assert( v.getdimension() == mesh.count_simplices(topdim) * (mesh.getinnerdimension()+1) );
    assert( v.isfinite() );

    auto sigmas = generateSigmas( IndexRange(1,1), IndexRange(0,topdim) );

    DenseMatrix conversion( topdim+1, topdim+1, 0 );
    
    for( int col = 0; col <= topdim; col++ ) {
        conversion( sigmas[col][1], col ) = 1.;
    }

    auto datafunction = [&](int c) -> FloatVector {
    
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        
        FloatVector extracted_coefficients( topdim+1 );
        
        for( int i = 0; i <= mesh.getinnerdimension(); i++ )
            extracted_coefficients[i] = v.at( c * (topdim+1) + i );
        
        FloatVector directions = mesh.getGradientMatrix(topdim,c) * ( conversion * extracted_coefficients );
        directions *= scaling; 

        return directions;
    };
    
    writeCellVectorData( datafunction, name, scaling );
    
    return *this;
}



VTKWriter VTKWriter::writeCellVectorData_barycentriccrosses( const FloatVector& v, const std::string& name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();

    assert( topdim == 3 );
    
    assert( v.getdimension() == 6 * mesh.count_simplices(topdim) );
    assert( v.isfinite() );

    auto sigmas = generateSigmas( IndexRange(1,2), IndexRange(0,topdim) );

    DenseMatrix conversion( 3, 6, 0 );
    
    for( int col = 0; col < 6; col++ ){
        
        // convert the spanning set coefficients to a basis
        // 0 1 -> 01 = - 21 - 31 =   12 + 13
        // 0 2 -> 02 = - 12 - 32 = - 12      + 23
        // 0 3 -> 03 = - 13 - 23 =      - 13 - 23 
        // 1 2                   =   12 
        // 1 3                   =        13
        // 2 3                   =             23
        FloatVector columnvalues(3);
               if( sigmas[col][1] == 0 and sigmas[col][2] == 1 ){
            columnvalues = FloatVector{  0.,  1.,  1. };            // +12+13
        } else if( sigmas[col][1] == 0 and sigmas[col][2] == 2 ){
            columnvalues = FloatVector{ -1.,  0.,  1. };            // -12+23
        } else if( sigmas[col][1] == 0 and sigmas[col][2] == 3 ){
            columnvalues = FloatVector{  0., -1., -1. };            // -13-23
        } else if( sigmas[col][1] == 1 and sigmas[col][2] == 2 ){
            columnvalues = FloatVector{  1.,  0.,  0. };
        } else if( sigmas[col][1] == 1 and sigmas[col][2] == 3 ){
            columnvalues = FloatVector{  0.,  1.,  0. };
        } else if( sigmas[col][1] == 2 and sigmas[col][2] == 3 ){
            columnvalues = FloatVector{  0.,  0.,  1. };
        } else {
            LOG << sigmas[col] << nl;
            unreachable();
        }
        conversion.setcolumn( col, columnvalues );

    }
    
    auto datafunction = [&](int c) -> FloatVector {
    
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        
        // extract the local coefficients 
        FloatVector extracted_coefficients( 6 );
        
        for( int i = 0; i < 6; i++ ) 
            extracted_coefficients[i] = v.at( c * 6 + i );
        
        // compute the crossproducts and the set the matrix 
        const auto J = mesh.getGradientMatrix(topdim,c);

        FloatVector cross_product_01 = FloatVector( { 
            J(1,0) * J(2,1) - J(2,0) * J(1,1), 
            J(2,0) * J(0,1) - J(0,0) * J(2,1), 
            J(0,0) * J(1,1) - J(1,0) * J(0,1)
        });
        
        FloatVector cross_product_02 = FloatVector( { 
            J(1,0) * J(2,2) - J(2,0) * J(1,2), 
            J(2,0) * J(0,2) - J(0,0) * J(2,2), 
            J(0,0) * J(1,2) - J(1,0) * J(0,2)
        });
        
        FloatVector cross_product_12 = FloatVector( { 
            J(1,1) * J(2,2) - J(2,1) * J(1,2), 
            J(2,1) * J(0,2) - J(0,1) * J(2,2), 
            J(0,1) * J(1,2) - J(1,1) * J(0,2)
        });

        DenseMatrix C(3,3);
        C.setcolumn( 0, cross_product_01 );
        C.setcolumn( 1, cross_product_02 );
        C.setcolumn( 2, cross_product_12 );
        
        // compute the direction 
        FloatVector directions = C * ( conversion * extracted_coefficients );
        directions *= scaling;

        return directions;
    };
    
    writeCellVectorData( datafunction, name, scaling );
    
    return *this;
}


VTKWriter VTKWriter::writeCellVectorData_Euclidean( int outerdim, const FloatVector& v, const std::string& name, Float scaling )
{
    const int topdim = mesh.getinnerdimension();
    
    assert( v.getdimension() == mesh.count_simplices(topdim) * (mesh.getinnerdimension()+1) );
    assert( v.isfinite() );
   
    auto datafunction = [&](int c) -> FloatVector {    
        
        assert( 0 <= c and c < this->mesh.count_simplices(topdim) );
        
        return v.getslice( outerdim*c, outerdim );
        
    };
    
    writeCellVectorData( datafunction, name, scaling );
    
    return *this;
}




VTKWriter VTKWriter::writePointCloud( const DenseMatrix& coords )
{
    assert( current_stage <= Stage::fielddata );
    current_stage = Stage::appendix;
    
    assert( coords.getdimin() == 2 or coords.getdimin() == 3 );
    assert( coords.isfinite() );


    os << "DATASET POLYDATA" << nl;
    os << "POINTS " << coords.getdimout() << " double" << nl;
    
    for( int v = 0; v < coords.getdimout(); v++ ) {
        os << space << coords(v,0);
        os << space << coords(v,1);
        if( coords.getdimin() == 3 ) 
            os << space << coords(v,2);
        else 
            os << space << 0.;
        os << nl;
    }
    
    os << nl;

    return *this;
}

    
