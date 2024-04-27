
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <functional>
#include <limits>
#include <memory>
#include <new>
#include <utility>
#include <vector>

#include "densematrix.hpp"

#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/simpleoperators.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../utility/random.hpp"


DenseMatrix::DenseMatrix( const DenseMatrix& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ), entries( new (std::nothrow) Float[ mat.getdimout() * mat.getdimin() ] ) // wierd
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = mat(r,c);
    DenseMatrix::check();
}
        
DenseMatrix::DenseMatrix( DenseMatrix&& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ), entries( std::move(mat.entries) )
{
    assert( entries != nullptr );
    mat.entries = nullptr;
    DenseMatrix::check();
}
        
DenseMatrix& DenseMatrix::operator=( const DenseMatrix& mat )
{
    assert( entries     != nullptr );
    assert( mat.entries != nullptr );
    assert( getdimin()  == mat.getdimin()  );
    assert( getdimout() == mat.getdimout() );
    
    if( this != &mat ) {
        for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            (*this)(r,c) = mat(r,c);
    }
    
    DenseMatrix::check();
    return *this;
}
        
DenseMatrix& DenseMatrix::operator=( DenseMatrix&& mat ) 
{
//     assert( entries     != nullptr );
    assert( mat.entries != nullptr );
    assert( getdimin()  == mat.getdimin()  );
    assert( getdimout() == mat.getdimout() );

    if( this != &mat ){
        delete[] entries;
        entries = mat.entries;
        mat.entries = nullptr;
    }

    DenseMatrix::check();
    return *this;
}
        
        
        
DenseMatrix::DenseMatrix( int dim, Float value )
: DenseMatrix( dim, dim, value )
{
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int dim, const std::function<Float(int,int)>& generator )
: DenseMatrix( dim, dim, generator )
{
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int dim, const std::vector<FloatVector>& coldata )
: DenseMatrix( dim, dim, coldata )
{
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int rows, int columns, Float value )
: LinearOperator( rows, columns ), entries( new (std::nothrow) Float[ rows * columns ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = value;
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int rows, int columns, const std::function<Float(int,int)>& generator )
: LinearOperator( rows, columns ), entries( new (std::nothrow) Float[ rows * columns ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = generator(r,c);
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int rows, int columns, const std::vector<FloatVector>& coldata )
: LinearOperator( rows, columns ), entries( new (std::nothrow) Float[ rows * columns ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = coldata.at(c).at(r);
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( int rows, int columns, const std::initializer_list<Float>& rowdata )
: LinearOperator( rows, columns ), entries( new (std::nothrow) Float[ rows * columns ] )
{
    assert( rowdata.size() == rows * columns );
    assert( entries != nullptr );

    int i = 0;
    for( const Float& f : rowdata ) entries[i++] = f;
    DenseMatrix::check();
}



DenseMatrix::DenseMatrix( const ScalingOperator& scaling )
: LinearOperator( scaling.getdimout(), scaling.getdimin() ), 
  entries( new (std::nothrow) Float[ scaling.getdimout() * scaling.getdimin() ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = ( ( r == c ) ? scaling.getscaling() : 0. );
    DenseMatrix::check();
}
        
DenseMatrix::DenseMatrix( const DiagonalOperator& dia )
: LinearOperator( dia.getdimout(), dia.getdimin() ), 
  entries( new (std::nothrow) Float[ dia.getdimout() * dia.getdimin() ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = ( ( r == c ) ? dia.getdiagonal().at(r) : 0. );
    DenseMatrix::check();
}
        
DenseMatrix::DenseMatrix( const SparseMatrix& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ), 
  entries( new (std::nothrow) Float[ matrix.getdimout() * matrix.getdimin() ]() )
{
    assert( entries != nullptr );
    for( const SparseMatrix::MatrixEntry& entry : matrix.getentries() )
    {
        (*this)( entry.row, entry.column ) += entry.value;
    }
    DenseMatrix::check();
}
        
DenseMatrix::DenseMatrix( const FloatVector& myvector )
: LinearOperator( myvector.getdimension(), 1 ), 
  entries( new (std::nothrow) Float[ myvector.getdimension() ] )
{
    assert( entries != nullptr );
    for( int r = 0; r < myvector.getdimension(); r++ )
    {
        (*this)( r, 0 ) = myvector[r];
    }
    DenseMatrix::check();
}


DenseMatrix::DenseMatrix( int number_of_blocks, const DenseMatrix& mat, Float scaling )
: LinearOperator( mat.getdimout() * number_of_blocks, mat.getdimin() * number_of_blocks), 
  entries( new (std::nothrow) Float[ number_of_blocks * number_of_blocks * mat.getdimout() * mat.getdimin() ] ) // wierd
{
    assert( number_of_blocks >= 0);
    assert( entries != nullptr );
    
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)( r, c ) = 0.;
    
    for( int b = 0; b < number_of_blocks; b++ )
    for( int r = 0; r < mat.getdimout(); r++ )
    for( int c = 0; c < mat.getdimin(); c++ )
        (*this)( b * mat.getdimout() + r, b * mat.getdimin() + c ) = scaling * mat( r, c );
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( const DenseMatrix& mat, Float scaling )
: LinearOperator( mat.getdimout(), mat.getdimin() ), entries( new (std::nothrow) Float[ mat.getdimout() * mat.getdimin() ] ) // wierd
{
    assert( entries != nullptr );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = scaling * mat(r,c);
    DenseMatrix::check();
}

DenseMatrix::DenseMatrix( DenseMatrix&& mat, Float scaling )
: LinearOperator( mat.getdimout(), mat.getdimin() ), entries( std::move(mat.entries) )
{
    assert( entries != nullptr );
    mat.entries = nullptr;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) *= scaling;
    DenseMatrix::check();
}

        
DenseMatrix::~DenseMatrix()
{
    // DenseMatrix::check(); // explicitly disabled, might be in moved-from state 
    
    if( entries != nullptr ){
        delete[] entries;
    }

}




void DenseMatrix::check() const 
{
    #ifdef NDEBUG
    return;
    #endif
    
    LinearOperator::check();
//     assert( entries != nullptr );
}

std::string DenseMatrix::text() const
{
    check();
    std::string ret = "Dense Matrix " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + nl;
    for( int r = 0; r < getdimout(); r++ ) {
        for( int c = 0; c < getdimin(); c++ )
            ret += printf_into_string( "% .10Le", (long double)at(r,c) ) + "\t";
        ret += nl;
    }
    return ret;
}

std::string DenseMatrix::data_as_text( bool indexed, bool print_as_list ) const
{
    const int nc_precision = 10;

    const int nc_width = 7 + nc_precision;
    
    std::string ret;

    if( print_as_list ){
        
        for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            if(indexed) 
                ret += printf_into_string("%*d %*d %*.*Le\n", nc_width, r, nc_width, c, nc_width, nc_precision, (long double)(*this)(r,c) );
            else
                ret += printf_into_string("%*.*Le\n", nc_width, nc_precision, (long double)(*this)(r,c) );
        
    } else {

        for( int r = 0; r < getdimout(); r++ ) {

            if(indexed) ret += printf_into_string( "%*d : ", nc_width, r );
            
            for( int c = 0; c < getdimin(); c++ )
                ret += printf_into_string( "%*.*Le ", nc_width, nc_precision, (long double)(*this)(r,c) );
            
            ret += '\n'; 
        }
            }

    return ret;
}


DenseMatrix DenseMatrix::clone() const
{
    return DenseMatrix( *this );
}




void DenseMatrix::apply( FloatVector& dest, const FloatVector& add, Float scaling ) const 
{
    add.check();
    dest.check();
    assert( dest.getdimension() == getdimout() );
    assert( add.getdimension() == getdimin() );
    
    dest.zero();
    
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        dest[r] += scaling * (*this)(r,c) * add[c];
    
    /*
    load r 0
    .run_outer_loop
    if r == dimout goto .end_outer_loop
      load c 0
      .run_inner_loop
      if c == dimin goto .end_inner_loop
        // perform multiplication
      inc
      goto .run_inner_loop
      .end_inner_loop
      inc r
      goto run_outer_loop
    .end_outer_loop
    */
    
    /*
    load r 0
    .run_outer_loop
        load c 0
        .run_inner_loop
        
        // perform with r and c
        
        inc c
        if c < dimin goto .run_inner_loop
    inc r
    if r < dimout goto .run_outer_loop
    */

}

Float DenseMatrix::get( int r, int c ) const
{
    return (*this)( r, c );
}

void DenseMatrix::set( int r, int c, Float v )
{
    (*this)( r, c ) = v;
}        

Float& DenseMatrix::at( int r, int c ) &
{
    assert( 0 <= r && r < getdimout() );
    assert( 0 <= c && c < getdimin() );
    assert( 0 <= r * getdimin() + c );
    assert( r * getdimin() + c < getdimin()*getdimout() );
    return entries[  r * getdimin() + c ];
}

const Float& DenseMatrix::at( int r, int c ) const &
{
    assert( 0 <= r && r < getdimout() );
    assert( 0 <= c && c < getdimin() );
    assert( 0 <= r * getdimin() + c );
    assert( r * getdimin() + c < getdimin()*getdimout() );
    return entries[  r * getdimin() + c ];
}

Float& DenseMatrix::operator()( int r, int c ) &
{
    assert( 0 <= r && r < getdimout() );
    assert( 0 <= c && c < getdimin() );
    assert( 0 <= r * getdimin() + c );
    assert( r * getdimin() + c < getdimin()*getdimout() );
    return entries[ r * getdimin() + c ];
}

const Float& DenseMatrix::operator()( int r, int c ) const &
{
    assert( 0 <= r && r < getdimout() );
    assert( 0 <= c && c < getdimin() );
    assert( 0 <= r * getdimin() + c );
    assert( r * getdimin() + c < getdimin()*getdimout() );
    return entries[ r * getdimin() + c ];
}

DenseMatrix DenseMatrix::submatrix( const IndexMap& rows, const IndexMap& columns ) const
{
    check();
    rows.check(); 
    columns.check();
    
    assert( rows.getSourceRange().min() == 0 );
    assert( rows.getSourceRange().max() <= getdimout() - 1 );
    assert( rows.getTargetRange().min() == 0 );
    assert( rows.getTargetRange().max() == getdimout() - 1 );
    assert( columns.getSourceRange().min() == 0 );
    assert( columns.getSourceRange().max() <= getdimin() - 1 );
    assert( columns.getTargetRange().min() == 0 );
    assert( columns.getTargetRange().max() == getdimin() - 1 );

    assert( rows.isstrictlyascending()    );
    assert( columns.isstrictlyascending() );
    
    DenseMatrix ret( rows.getSourceRange().cardinality(), columns.getSourceRange().cardinality(), notanumber );
    
    for( int nr = 0; nr < ret.getdimout(); nr++ )
    for( int nc = 0; nc < ret.getdimin();  nc++ )
    {
        assert( rows.getSourceRange().contains(nr) );
        assert( columns.getSourceRange().contains(nc) );
        
        ret( nr, nc ) = (*this)( rows[nr], columns[nc] );
    }
    
    return ret;
    
//     if( getdimin() == 0 || getdimout() == 0 )
//         return *this;
//     
//     assert( rows.getSourceRange().min() == 0 );
//     assert( rows.getSourceRange().max() <= getdimout() - 1 );
//     assert( rows.getTargetRange().min() == 0 );
//     assert( rows.getTargetRange().max() == getdimout() - 1 );
//     assert( columns.getSourceRange().min() == 0 );
//     assert( columns.getSourceRange().max() <= getdimin() - 1 );
//     assert( columns.getTargetRange().min() == 0 );
//     assert( columns.getTargetRange().max() == getdimin() - 1 );
//     assert( rows.isstrictlyascending() );
//     assert( columns.isstrictlyascending() );
//     
//     DenseMatrix ret( maximum(0,rows.getSourceRange().max() - 1), 
//                     maximum(0,columns.getSourceRange().max() - 1) );
//     for( int nr = 0; nr < ret.getdimout(); nr++ )
//     for( int nc = 0; nc < ret.getdimin();  nc++ )
//         ret( nr, nc ) = (*this)( rows[nr], columns[nc] );
//     
//     ret.check();
//     return ret;
    
}

FloatVector DenseMatrix::getDiagonal() const {
    int dim = minimum( getdimin(), getdimout() );
    FloatVector ret( dim );
    for( int d = 0; d < dim; d++ ) ret[d] = at(d,d);
    return ret;
}
        
        


FloatVector DenseMatrix::flattencolumns() const
{
    check();
    FloatVector ret( getdimin() * getdimout() );
    for( int c = 0; c < getdimin();  c++ )
    for( int r = 0; r < getdimout(); r++ )
        ret[ c * getdimout() + r ] = get(r,c);
    return ret;
}

FloatVector DenseMatrix::flattenrows() const
{
    check();
    FloatVector ret( getdimin() * getdimout() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        ret[ r * getdimin() + c ] = get(r,c);
    return ret;
}

// std::vector<SparseMatrix::MatrixEntry> DenseMatrix::getSparseMatrixEntries( bool clean ) const
// {
//     check();
//     std::vector<SparseMatrix::MatrixEntry> ret(0);
//     ret.reserve( getdimin() * getdimout() );
//     for( int r = 0; r < getdimout(); r++ )
//     for( int c = 0; c < getdimin (); c++ )
//         if( not clean or get(r,c) != 0. )
//             ret.push_back( { r, c, get(r,c) } );
//     ret.reserve( ret.size() );
//     return ret;
// }









void DenseMatrix::zeromatrix()
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = 0.;
}

void DenseMatrix::unitmatrix()
{
    check();
    assert( getdimout() == getdimin() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        if( r == c ) 
            (*this)(r,c) = 1.;
        else
            (*this)(r,c) = 0.;
}

void DenseMatrix::randommatrix()
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = gaussrand();
}

void DenseMatrix::randomintegermatrix( int min, int max )
{
    check();
    assert( min < max );
    assert( min+2 <= max );
    unsigned long int diff = max - min;
    assert( diff+1 > 0 );
    LOG << diff << nl << max << nl << min << nl;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
    {
        unsigned int random_int = random_integer();
        Float value = random_int % (diff+1) + (Float)min; // force the addition to be floating-point, avoid unsigned arithmetic
        (*this)(r,c) = value;
    }
        
}

void DenseMatrix::random_orthogonal_matrix()
{
    check();
    assert( issquare() );

    const auto N = getdimin();
    auto rows = std::vector<FloatVector>( N, FloatVector(N) );
    
    for( auto& row : rows ) { 
        row.random(); 
        row.normalize(); 
        assert( row.isfinite() ); 
        assert( std::isfinite( row.norm() ) ); 
    }
    
    for( int i = 0; i < N; i++ ) {
        for( int j = 0; j < i; j++ )
            rows[i] -= (rows[i]*rows[j]) * rows[j];
        rows[i].normalize();
    }
    
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = rows[r][c];
}

void DenseMatrix::scale( Float s )
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) *= s;
}

void DenseMatrix::set( Float s )
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) = s;
}

void DenseMatrix::add( Float s )
{
    check();
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r,c) += s;
}

FloatVector DenseMatrix::getrow( int r ) const 
{
    check();
    assert( 0 <= r && r < getdimout() );
    FloatVector row( getdimin() );
    for( int c = 0; c < getdimin(); c++ )
            row[c] = (*this)(r,c);
    row.check();
    return row;
}

void DenseMatrix::setrow( int r, const FloatVector& row )
{
    check();
    row.check();
    assert( 0 <= r && r < getdimout() );
    assert( row.getdimension() == getdimin() );
    for( int c = 0; c < getdimin(); c++ )
      (*this)(r,c) = row[c];
}

FloatVector DenseMatrix::getcolumn( int c ) const 
{
    check();
    assert( 0 <= c && c < getdimin() );
    FloatVector column( getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        column[r] = (*this)(r,c);
    column.check();
    return column;
}

void DenseMatrix::setcolumn( int c, const FloatVector& column )
{
    check();
    column.check();
    assert( 0 <= c && c < getdimin() );
    assert( column.getdimension() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        (*this)(r,c) = column[r];
}

void DenseMatrix::addrow( int r, const FloatVector& row, Float s )
{
    check();
    row.check();
    assert( 0 <= r && r < getdimout() );
    assert( row.getdimension() == getdimin() );
    for( int c = 0; c < getdimin(); c++ )
      (*this)(r,c) += s * row[c];
}

void DenseMatrix::addcolumn( int c, const FloatVector& column, Float s )
{
    check();
    column.check();
    assert( 0 <= c && c < getdimin() );
    assert( column.getdimension() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        (*this)(r,c) += s * column[r];
}

void DenseMatrix::swaprow( int r1, int r2 )
{
    check();
    assert( 0 <= r1 && r1 < getdimout() );
    assert( 0 <= r2 && r2 < getdimout() );
    if( r1 == r2 ) return;
    for( int c = 0; c < getdimin(); c++ )
        std::swap( (*this)(r1,c), (*this)(r2,c) );
    check();
}

void DenseMatrix::swapcolumn( int c1, int c2 )
{
    check();
    assert( 0 <= c1 && c1 < getdimin() );
    assert( 0 <= c2 && c2 < getdimin() );

    if( c1 == c2 ) return;

    for( int r = 0; r < getdimout(); r++ )
        std::swap( (*this)(r,c1), (*this)(r,c2) );

    check();
}


void DenseMatrix::scalerow( int r, Float alpha )
{
    check();
    assert( 0 <= r && r < getdimout() );
    for( int c = 0; c < getdimin(); c++ )
      (*this)(r,c) *= alpha;
}

void DenseMatrix::scalecolumn( int c, Float alpha )
{
    check();
    assert( 0 <= c && c < getdimin() );
    for( int r = 0; r < getdimout(); c++ )
      (*this)(r,c) *= alpha;
}

void DenseMatrix::addrow( int r1, int r2, Float alpha )
{
    check();
    assert( 0 <= r1 && r1 < getdimout() );
    assert( 0 <= r2 && r2 < getdimout() );
    for( int c = 0; c < getdimin(); c++ )
        (*this)(r1,c) += alpha * (*this)(r2,c);
    check();
}

void DenseMatrix::addcolumn( int c1, int c2, Float alpha )
{
    check();
    assert( 0 <= c1 && c1 < getdimin() );
    assert( 0 <= c2 && c2 < getdimin() );
    for( int r = 0; r < getdimout(); r++ )
        (*this)(r,c1) += alpha * (*this)(r,c2);
    check();
}


        




        
void DenseMatrix::indexmapping( const IndexMap& im )
{
    check();
    im.check();
    if( im.getSourceRange().isempty() )
        return;
    assert( im.getSourceRange().min() == 0 );
    assert( im.getSourceRange().max() == getdimin() - 1 );
    assert( im.getTargetRange().min() == 0 );
    assert( im.getTargetRange().max() == getdimout() - 1 );
    zeromatrix();
    for( int c = 0; c < im.getSourceRange().max() - 1; c++ )
        (*this)( im[c], c ) = 1.;
    check();
}





        
void DenseMatrix::add( const DenseMatrix& addendum )
{
    check();
    addendum.check();
    for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            (*this)(r,c) = (*this)(r,c) + addendum(r,c);
    check();
}

void DenseMatrix::add( Float t, const DenseMatrix& addendum )
{
    check();
    addendum.check();
    for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            (*this)(r,c) = (*this)(r,c) + t * addendum(r,c);
    check();
}

void DenseMatrix::add( Float s, Float t, const DenseMatrix& addendum )
{
    check();
    addendum.check();
    for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            (*this)(r,c) = s * (*this)(r,c) + t * addendum(r,c);
    check();
}
        



DenseMatrix DenseMatrix::symmetricPart() const 
{
    assert( issquare() );
    DenseMatrix ret( getdimout() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret(r,c) = 0.5 * (*this)(r,c) + 0.5 * (*this)(c,r);
    return ret;
}

DenseMatrix DenseMatrix::antisymmetricPart() const 
{
    assert( issquare() );
    DenseMatrix ret( getdimout() );
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret(r,c) = 0.5 * (*this)(r,c) - 0.5 * (*this)(c,r);
    return ret;
}






Float DenseMatrix::maxabsoluteentry() const
{
    check();
    Float ret = - std::numeric_limits<Float>::infinity();
    for( int r = 0; r < getdimout(); r++ )
        for( int c = 0; c < getdimin(); c++ )
            ret = maximum( ret, (*this)(r,c) );
    return ret;
}









bool DenseMatrix::issquare() const
{
    return getdimout() == getdimin();
}

bool DenseMatrix::issymmetric() const
{
    if( not issquare() ) return false;
    
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c <= r; c++ )
        if( not is_numerically_close( get(r,c), get(c,r) ) ) 
            return false;
    
    return true;
}

bool DenseMatrix::isantisymmetric() const
{
    if( not issquare() ) return false;
    
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c <= r; c++ )
        if( not is_numerically_close( get(r,c), -get(c,r) ) ) 
            return false;
    
    return true;
}

bool DenseMatrix::isfinite() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( not std::isfinite( get(r,c) ) ) 
            return false;
    return true;
}

bool DenseMatrix::iszero() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) != 0. ) 
            return false;
    return true;
}



bool DenseMatrix::ispositive() const 
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) <= 0. ) 
            return false;
    return true;
}

bool DenseMatrix::isnegative() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) >= 0. ) 
            return false;
    return true;
}

bool DenseMatrix::isnonnegative() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) < 0. ) 
            return false;
    return true;
}

bool DenseMatrix::isnonpositive() const
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        if( get(r,c) > 0. ) 
            return false;
    return true;
}

bool DenseMatrix::is_numerically_small( Float threshold ) const
{
    return this->sumnorm() < threshold;
}

bool DenseMatrix::is_numerically_identity( Float threshold ) const
{
    //Float sum = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
    {
        if ( r == c and absolute( get(r,c) - 1.0 ) >= threshold ) return false;

        if ( r != c and absolute( get(r,c) ) >= threshold ) return false;
    }
    return true;
}


Float DenseMatrix::norm() const
{
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret += get(r,c) * get(r,c);
    return std::sqrt( ret );
}

Float DenseMatrix::frobeniusnorm() const
{
    return norm();
}

Float DenseMatrix::maxnorm() const
{
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret = maximum( ret, absolute( get(r,c) ) );
    return ret;
}

Float DenseMatrix::sumnorm() const
{
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret += absolute( get(r,c) );
    return ret;
}

Float DenseMatrix::lpnorm( Float p ) const
{
    assert( p >= 1. );
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin();  c++ )
        ret += power_numerical( absolute( get(r,c) ), p );
    return power_numerical( ret, 1. / p );
}

Float DenseMatrix::norm_row_col( Float p, Float q ) const 
{
    check();
    assert( 1. <= p && 1. <= q );
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ ) {
        Float zeile = 0.;
        for( int c = 0; c < getdimin(); c++ )
            zeile += power_numerical( absolute( get(r,c) ), q );
        ret += power_numerical( zeile, p/q );
    }
    ret = power_numerical( ret, 1. / p );
    return ret;
}

Float DenseMatrix::norm_col_row( Float p, Float q ) const 
{
    check();
    assert( 1. <= p && 1. <= q );
    Float ret = 0.;
    for( int c = 0; c < getdimin(); c++ ) {
        Float spalte = 0.;
        for( int r = 0; r < getdimout(); r++ )
            spalte += power_numerical( absolute( get(r,c) ), q );
        ret += power_numerical( spalte, p/q );
    }
    ret = power_numerical( ret, 1. / p );
    return ret;
}

Float DenseMatrix::NormOperatorL1() const 
{
    check();
    Float ret = 0.;
    for( int c = 0; c < getdimin(); c++ ) {
        Float spalte = 0.;
        for( int r = 0; r < getdimout(); r++ )
            spalte += absolute( get(r,c) );
        ret = maximum( ret, spalte );
    }
    return ret;
}

Float DenseMatrix::NormOperatorMax() const 
{
    check();
    Float ret = 0.;
    for( int r = 0; r < getdimout(); r++ ) {
        Float zeile = 0.;
        for( int c = 0; c < getdimin(); c++ )
            zeile += absolute( get(r,c) );
        ret = maximum( ret, zeile );
    }
    return ret;
}






Float DenseMatrix::trace() const 
{
    check();
    assert( issquare() );
    Float ret = 0.;
    for( int i = 0; i < getdimout(); i++ )
      ret += get(i,i);
    return ret;
}


DenseMatrix DenseMatrix::Gerschgorin() const 
{
    check();
    assert( issquare() );
    return GerschgorinRow();
}

DenseMatrix DenseMatrix::GerschgorinRow() const 
{
    check();
    assert( issquare() );
    DenseMatrix ret( getdimout(), 2 );
    for( int r = 0; r < getdimout(); r++ )
    {
        ret( r, 0 ) = get(r,r);
        ret( r, 1 ) = 0.;
        for( int c = 0; c < r; c++ )
            ret( r, 1 ) += absolute( get(r,c) );
        for( int c = r+1; c < getdimout(); c++ )
            ret( r, 1 ) += absolute( get(r,c) );
    }
    return ret;
}

DenseMatrix DenseMatrix::GerschgorinColumn() const
{
    check();
    assert( issquare() );
    DenseMatrix ret( getdimout(), 2 );
    for( int c = 0; c < getdimout(); c++ )
    {
        ret( c, 0 ) = get(c,c);
        ret( c, 1 ) = 0.;
        for( int r = 0; r < c; r++ )
            ret( c, 1 ) += absolute( get(r,c) );
        for( int r = c+1; r < getdimout(); r++ )
            ret( c, 1 ) += absolute( get(r,c) );
    }
    return ret;
}

Float DenseMatrix::EigenvalueEstimate() const 
{
    assert( issquare() );
    return maxnorm() * getdimout();
}







Float* DenseMatrix::raw()
{
    return entries;
}

const Float* DenseMatrix::raw() const
{
    return entries;
}

/* Memory size */
        
std::size_t DenseMatrix::memorysize() const
{
    return sizeof(*this) + getdimin() * getdimout() * sizeof(this->entries[0]);
}

















DenseMatrix IdentityMatrix( int dim )
{
    return DenseMatrix( dim, []( int r, int c ) -> Float{ return r==c ? 1. : 0.; } );
}




DenseMatrix MatrixMult( const DenseMatrix& left, const DenseMatrix& right )
{
    left.check();
    right.check();

    const int lin = left.getdimin();
    const int lout = left.getdimout();
    const int rin = right.getdimin();
    const int rout = right.getdimout();

    assert( lin == rout );

    DenseMatrix ret( lout, rin, 0. );
    
    for( int lo = 0; lo < lout; lo++ )
    for( int ri = 0; ri < rin; ri++ )
    for( int m = 0; m < rout; m++ )
        ret( lo, ri ) += left( lo, m ) * right( m, ri );

    ret.check();
    return ret;
}


DenseMatrix MatrixTripleMult( const DenseMatrix& A, const DenseMatrix& B )
{
    A.check();
    B.check();

    assert( A.issquare());

    const int a = A.getdimin();
    const int b = B.getdimin();
    
    assert( a == B.getdimout() );

    DenseMatrix ret( b, b, 0. );

    for( int j = 0; j < b; j++ )
    for( int k = 0; k < a; k++ ) 
    {
        Float AB_kj = 0;
    
        for( int l = 0; l < a; l++ )
            AB_kj += A(k,l) * B(l,j);
    
        for( int i = 0; i < b; i++ )
            ret(i,j) += B(k,i) * AB_kj;
    }


    // for( int i = 0; i < b; i++ )
    // for( int k = 0; k < a; k++ )
    // for( int j = 0; j < b; j++ )
    // for( int l = 0; l < a; l++ )
    //     ret( i, j ) += B(k,i) * A(k,l) * B(l,j);

    ret.check();
    return ret;
}


DenseMatrix HilbertMatrix( int n )
{
    std::function<Float(int,int)> hilbertmatrix_generator = [](int r, int c) { return 1. / (r+c+1); };

    return DenseMatrix( n, hilbertmatrix_generator );
}

DenseMatrix InvHilbertMatrix( int n )
{
    std::function<Float(int,int)> invhilbertmatrix_generator = [=](int r, int c) { 
        // https://mathoverflow.net/questions/47561/deriving-inverse-of-hilbert-matrix
        const int i = r+1;
        const int j = c+1;
        return sign_power(i+j) * (i+j-1)
               * binomial_integer( n+i-1, n-j )
               * binomial_integer( n+j-1, n-i )
               * square( binomial_integer(i+j-2,i-1) );
    };

    return DenseMatrix( n, invhilbertmatrix_generator );
}
