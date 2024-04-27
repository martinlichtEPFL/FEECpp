
#include <cstdlib>

#include <algorithm>
#include <functional>
#include <initializer_list>
#include <ostream>
#include <utility>
#include <vector>

#include "sparsematrix.hpp"
#include "../dense/densematrix.hpp"

const bool sparse_matrix_verbosity = false;


SparseMatrix::SparseMatrix( int dimout, int dimin, int numentries, std::function<MatrixEntry(int)> generator )
: LinearOperator( dimout, dimin ), 
  entries(0)
{
    assert( 0 <= numentries );
    assert( numentries <= entries.max_size() );
    
    if( sparse_matrix_verbosity ) LOG << "SparseMatrix allocation: " << numentries << nl;
    
    entries.reserve( numentries );
    
    for( int i = 0; i < numentries; i++ ) 
    {
      SparseMatrix::MatrixEntry entry = generator(i);
      
      assert( 0 <= entry.row && entry.row < getdimout() );
      assert( 0 <= entry.column && entry.column < getdimin() );
      
      entries.push_back( entry );
    }
    
    SparseMatrix::check();    
}

// SparseMatrix::SparseMatrix( int dimout, int dimin )
// : SparseMatrix( dimout, dimin, std::vector<MatrixEntry>(0) )
// {
//     SparseMatrix::check();
// }

SparseMatrix::SparseMatrix( int dimout, int dimin, const std::vector<MatrixEntry>& entries )
: LinearOperator( dimout, dimin ), entries(entries)
{
    SparseMatrix::check();
}

SparseMatrix::SparseMatrix( int dimout, int dimin, const std::initializer_list<MatrixEntry>& ent )
: LinearOperator( dimout, dimin ), entries(0)
{
    entries.reserve( ent.size() );
    std::copy( ent.begin(), ent.end(), entries.begin() );
    assert( ent.size() == entries.size() );
}

SparseMatrix::SparseMatrix( const FloatVector& diagonal )
: LinearOperator( diagonal.getdimension() ),
  entries( diagonal.getdimension() )
{
    assert( getdimin() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        entries[r] = { r, r, diagonal.at(r) };
    SparseMatrix::check();    
}

SparseMatrix::SparseMatrix( const ScalingOperator& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ),
  entries( matrix.getdimout() )
{
    assert( getdimin() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        entries[r] = { r, r, matrix.getscaling() };
    SparseMatrix::check();    
}

SparseMatrix::SparseMatrix( const DiagonalOperator& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ),
  entries( matrix.getdimout() )
{
    assert( getdimin() == getdimout() );
    for( int r = 0; r < getdimout(); r++ )
        entries[r] = { r, r, matrix.getdiagonal().at(r) };
    SparseMatrix::check();    
}

SparseMatrix::SparseMatrix( const DenseMatrix& matrix )
: LinearOperator( matrix.getdimout(), matrix.getdimin() ), 
  entries( matrix.getdimout() )
{
    for( int r = 0; r < getdimout(); r++ )
    for( int c = 0; c < getdimin(); c++ )
        entries[ r * getdimin() + c ] = { r, c, matrix.get(r,c) };
    SparseMatrix::check();    
}

SparseMatrix::~SparseMatrix()
{
    SparseMatrix::check();    
}






SparseMatrix::SparseMatrix( const SparseMatrix& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ),
  entries( mat.entries )
{
    LOG << "*************************************************\n";
    LOG << "*********** WARNING: DEEP COPY ******************\n";
    LOG << "***********  OF SPARSE MATRIX  ******************\n";
    LOG << "*************************************************\n";
    SparseMatrix::check();
}

SparseMatrix& SparseMatrix::operator=( const SparseMatrix& mat )
{
    LOG << "*************************************************\n";
    LOG << "********** WARNING: DEEP ASSIGN *****************\n";
    LOG << "***********  OF SPARSE MATRIX  ******************\n";
    LOG << "*************************************************\n";
    assert( getdimin() == mat.getdimin() );
    assert( getdimout() == mat.getdimout() );
    this->entries = mat.entries;
    check();
    return *this;
}

SparseMatrix::SparseMatrix( SparseMatrix&& mat )
: LinearOperator( mat.getdimout(), mat.getdimin() ),
  entries( std::move(mat.entries) )
{
    SparseMatrix::check();
}

SparseMatrix& SparseMatrix::operator=( SparseMatrix&& mat )
{
    assert( getdimin() == mat.getdimin() );
    assert( getdimout() == mat.getdimout() );
    this->entries = std::move( mat.entries );
    check();
    return *this;
}











void SparseMatrix::check() const
{
    #ifdef NDEBUG
    return;
    #endif
    
    LinearOperator::check();
    for( const MatrixEntry& me : entries ) {
        assert( 0 <= me.row && me.row < getdimout() );
        assert( 0 <= me.column && me.column < getdimin() );
    }
}

std::string SparseMatrix::text() const
{
    std::string ret = "SparseMatrix " + std::to_string(getdimout()) + "x" + std::to_string(getdimin());
    for( const MatrixEntry& entry : entries )
        ret += ( "\n" + std::to_string(entry.row) + " " + std::to_string(entry.column) + " : " + std::to_string(entry.value) );
        // ret += ( "\n" + std::to_string(entry.row) + " " + std::to_string(entry.column) + " : " + printf_into_string("%.17le",entry.value) );
    return ret;
}

// void SparseMatrix::printplain( std::ostream& os ) const
// {
//     for( const MatrixEntry& entry : entries )
//         os << entry.row << " " << entry.column << " " << entry.value << nl;
//     os << nl;
// }


void SparseMatrix::apply( FloatVector& dest, const FloatVector& add, Float scaling ) const 
{
    check();
    add.check();
    dest.check();
    
    assert( getdimin() == add.getdimension() );
    assert( getdimout() == dest.getdimension() );

    dest.zero();
    
    // for( const MatrixEntry& rcv : entries )
    //     dest.setentry( rcv.row, dest.getentry( rcv.row ) + scaling * rcv.value * add.getentry( rcv.column ) );

    for( int e = 0; e < entries.size(); e++ ){
        
        dest[ entries[e].row ] += scaling * entries[e].value * add[ entries[e].column ];
        
    }


}










void SparseMatrix::scale ( Float s )
{
    for( auto& e : this->entries ) e.value *= s;
}

bool SparseMatrix::isfinite() const
{
    for( const MatrixEntry& rcv : entries )
        if( not std::isfinite( rcv.value ) )
            return false;
    return true;
}

FloatVector SparseMatrix::getDiagonal() const
{ 
    check();
    assert( getdimin() == getdimout() );
    auto ret = FloatVector( getdimin(), 0. );

    for( const auto& entry : entries )
        if( entry.row == entry.column )
            ret[ entry.row ] += entry.value;
    
    return ret;
}


int SparseMatrix::getnumberofzeroentries() const
{
    check();
    
    int ret = 0;
    
    for( const auto& entry : entries )
        if( entry.value == 0. )
            ret++;
    
    return ret;
    
}




        
        







const std::vector<SparseMatrix::MatrixEntry>& SparseMatrix::getentries() const
{
    return entries;
}

std::vector<SparseMatrix::MatrixEntry>& SparseMatrix::getentries()
{
    return entries;
}
        
int SparseMatrix::getnumberofentries() const 
{
    return SIZECAST( entries.size() );
}










bool SparseMatrix::is_sorted( SparseMatrix::MatrixEntrySorting manner ) const
{
    for( int i = 1; i < entries.size(); i++ ) {

        const auto& a = entries[i-1];
        const auto& b = entries[  i];

        if( manner == MatrixEntrySorting::rowwise ){
            if( a.row    > b.row    or ( a.row == b.row and a.column > b.column ) )
                return false;
        }else{
            if( a.column > b.column or ( a.column == b.column and a.row > b.row ) )
                return false;
        }

    }

    return true;
    
}


UNUSED static int internal_compare_rowfirst( const void* _a, const void* _b )
{
    const SparseMatrix::MatrixEntry* a = static_cast<const SparseMatrix::MatrixEntry*>(_a);
    const SparseMatrix::MatrixEntry* b = static_cast<const SparseMatrix::MatrixEntry*>(_b);
    
    if( a->row < b->row ) {
        return -1;
    } else if( a->row > b->row ) {
        return 1;
    } else if( a->column < b->column ) {
        return -1;
    } else if( a->column > b->column ) {
        return 1;
    } else {
        return 0;
    }
}
    
UNUSED static int internal_compare_colfirst( const void* _a, const void* _b )
{
    const SparseMatrix::MatrixEntry* a = static_cast<const SparseMatrix::MatrixEntry*>(_a);
    const SparseMatrix::MatrixEntry* b = static_cast<const SparseMatrix::MatrixEntry*>(_b);
    
    if( a->column < b->column ) {
        return -1;
    } else if( a->column > b->column ) {
        return 1;
    } else if( a->row < b->row ) {
        return -1;
    } else if( a->row > b->row ) {
        return 1;
    } else {
        return 0;
    }
}


const SparseMatrix& SparseMatrix::sortentries( SparseMatrix::MatrixEntrySorting manner ) const
{
    check();

    // if( manner == MatrixEntrySorting::rowwise )
    //     std::qsort( entries.data(), entries.size(), sizeof(MatrixEntry), &internal_compare_rowfirst );
    // else 
    //     std::qsort( entries.data(), entries.size(), sizeof(MatrixEntry), &internal_compare_colfirst );

    if( manner == MatrixEntrySorting::rowwise )
        std::sort( entries.begin(), entries.end(), []( const MatrixEntry& a, const MatrixEntry& b) {
            return a.row < b.row || ( a.row == b.row && a.column < b.column );   
        });
    else 
        std::sort( entries.begin(), entries.end(), []( const MatrixEntry& a, const MatrixEntry& b) {
            return a.column < b.column || ( a.column == b.column && a.row < b.row );   
        });
    
    
//     const int N = getnumberofentries();
//     
//     for( int i = 0; i < N; i++ ) {
//     
//         MatrixEntry entry = entries.at(i);
//         int j = i;
//         while( j > 0 and ( entries.at(j-1).row > entry.row or entries.at(j-1).column > entry.column ) ) {
//             entries.at(j) = entries.at(j-1);
//             j--;
//         }
//         entries.at(j) = entry;
//     }
    
    assert( is_sorted( manner ) );

    check();
    
    return *this;
}

const SparseMatrix& SparseMatrix::sortandcompressentries( SparseMatrix::MatrixEntrySorting manner ) const
{
    check();

    if( sparse_matrix_verbosity ) LOG << "SparseMatrix: Remove zeroes..." << nl;
    
    {
        
        // for( int c = 0; c < entries.size(); c++ ) 
        //     if( entries[c].value == 0 ) { entries.erase( entries.begin()+c ); c--; }
        
        auto it = std::remove_if(entries.begin(), entries.end(), 
                    []( const SparseMatrix::MatrixEntry& e ){ return e.value == 0.; } 
                  );
        entries.erase( it, entries.end() );
        
        for( int c = 0; c < entries.size(); c++ ) assert( entries[c].value != 0. );
        
    }
    
    if( sparse_matrix_verbosity ) LOG << "SparseMatrix: Sorting..." << nl;
    
    sortentries( manner );
    
    if( sparse_matrix_verbosity ) LOG << "SparseMatrix: Compressing..." << nl;
    
    assert( is_sorted(manner) );
    
    {
        
        int dest = 0;
        
        for( int src = 1; src < getnumberofentries(); src++ )
        {
            assert( dest < src );
        
            if( entries[src].row == entries[dest].row and entries[src].column == entries[dest].column ) {
            
                entries[dest].value += entries[src].value;
            
            } else {
            
                dest++;
                entries[dest] = entries[src];
            
            }
            
            assert( dest <= src );
        }
    
        entries.resize( dest+1 );
    
        check(); 
        
    }
    
    assert( is_sorted(manner) );
    
    if( sparse_matrix_verbosity ) LOG << "SparseMatrix: Done!" << nl;
    
    return *this;
}






        
        
        
        
void SparseMatrix::reserve( int n ) const
{
    entries.reserve( n );
}
        
        


const SparseMatrix::MatrixEntry& SparseMatrix::getentry( int i ) const
{
    assert( 0 <= i && i < getnumberofentries() );
    return entries[i];
}

SparseMatrix::MatrixEntry& SparseMatrix::getentry( int i )
{
    assert( 0 <= i && i < getnumberofentries() );
    return entries[i];
}
        
void SparseMatrix::setentry( int i, int r , int c, Float v )
{
//     check();
    SparseMatrix::MatrixEntry temp;
    temp.row = r;
    temp.column = c;
    temp.value = v;
    setentry( i, temp );
//     check();
}

void SparseMatrix::setentry( int i, MatrixEntry entry )
{
//     check();
    assert( 0 <= entry.row && entry.row <= getdimout() );
    assert( 0 <= entry.column && entry.column <= getdimin() );
    assert( 0 <= i && i < entries.size() );
    entries.at(i) = entry;
//     check();
}
        
        
void SparseMatrix::appendentry( int r, int c, Float v )
{
//     check();
    SparseMatrix::MatrixEntry temp;
    temp.row = r;
    temp.column = c;
    temp.value = v;
    appendentry( temp );
//     check();
}

void SparseMatrix::appendentry( SparseMatrix::MatrixEntry entry )
{
//     check();
    assert( 0 <= entry.row && entry.row <= getdimout() );
    assert( 0 <= entry.column && entry.column <= getdimin() );
    entries.push_back( entry );
//     check();
}

void SparseMatrix::clearentries()
{
    check();
    entries.clear();
    check();
}









SparseMatrix SparseMatrix::getTranspose() const 
{
    
    auto newentries = entries;

    for( auto& newentry : newentries )
        std::swap( newentry.row, newentry.column );

    return SparseMatrix( getdimin(), getdimout(), newentries );

    // SparseMatrix ret = SparseMatrix( getdimout(), getdimin() );
    
    // for( const MatrixEntry& me : entries ) {
    //     ret.entries.push_back( { me.column, me.row, me.value } );
    // }
    
    // return ret;
}



// FloatVector InverseDiagonalPreconditioner( const SparseMatrix& mat ) const
// { 
//     check();
//     assert( mat.getdimin() == mat.getdimout() );
//     auto ret = FloatVector( mat.getdimin(), 0. );
//     for( auto& entry : mat.entries )
//         if( entry.row == entry.column )
//             ret[ entry.row ] += entry.value;
//     for( int r = 0; r < mat.getdimout(); r++ ) {
//         assert( ret[r] >= 0. );
//         if( ret[r] > 0. ) ret[r] = 1. / ret[r];
//     }
//     return ret;
// }




/* Memory size */
        
std::size_t SparseMatrix::memorysize() const
{
    return sizeof(*this) + entries.size() * sizeof(decltype(entries)::value_type);
}










SparseMatrix SparseMatrixMultiplication( const SparseMatrix& left, const SparseMatrix& right )
{

    assert( left.getdimin() == right.getdimout() );
    
//     LOG << "--- SparseMatrix Product" << nl;
//     LOG << "--- Sort and compress" << nl;

    // SectionProfiler beacon;
    
    left.sortandcompressentries( SparseMatrix::MatrixEntrySorting::columnwise );
    right.sortandcompressentries( SparseMatrix::MatrixEntrySorting::rowwise );

    // beacon.ping("Sorted"); // LOG << "--- Counting" << nl;
    
    const int lnum = left.getnumberofentries();
    const int rnum = right.getnumberofentries();
    const SparseMatrix::MatrixEntry* const ldata =  left.getentries().data();
    const SparseMatrix::MatrixEntry* const rdata = right.getentries().data();
    int counter = 0;
    
    {

        int lbase = 0; int rbase = 0;
        
        while( lbase < lnum and rbase < rnum ){

            if( ldata[lbase].column == rdata[rbase].row ) {
                
                int index = ldata[lbase].column;
                assert( index == rdata[rbase].row );
                int lnext = lbase;
                int rnext = rbase;
                while( lnext < lnum and ldata[lnext].column == index ) lnext++;
                while( rnext < rnum and rdata[rnext].row    == index  ) rnext++;
                counter = counter + ( lnext - lbase ) * ( rnext - rbase );
                lbase = lnext;
                rbase = rnext;
                
            } else if( ldata[lbase].column < rdata[rbase].row ) {
                
                lbase++;
                
            } else if( rdata[rbase].row < ldata[lbase].column ) {
                
                rbase++;
                
            } else {
                unreachable();
            }

        }

    }
    
    // beacon.ping("Counted"); // LOG << "--- Assemble" << nl;

    std::vector<SparseMatrix::MatrixEntry> new_entries;
    new_entries.reserve( counter );
    
    {

        int lbase = 0; int rbase = 0;
        
        while( lbase < lnum and rbase < rnum ){

            if( ldata[lbase].column == rdata[rbase].row ) {
                
                int index = ldata[lbase].column;
                assert( index == rdata[rbase].row );
                int lnext = lbase;
                int rnext = rbase;
                while( lnext < lnum and ldata[lnext].column == index ) lnext++;
                while( rnext < rnum and rdata[rnext].row    == index ) rnext++;
                
                for( int lcurr = lbase; lcurr < lnext; lcurr++ )
                for( int rcurr = rbase; rcurr < rnext; rcurr++ )
                    new_entries.push_back( { 
                        ldata[lcurr].row, 
                        rdata[rcurr].column, 
                        ldata[lcurr].value * rdata[rcurr].value
                    } );
                
                lbase = lnext;
                rbase = rnext;
                
            } else if( ldata[lbase].column < rdata[rbase].row ) {
                
                lbase++;
                
            } else if( rdata[rbase].row < ldata[lbase].column ) {
                
                rbase++;
                
            } else {
                unreachable();
            }

        }

    }

    // beacon.ping("Assembled"); // LOG << "--- Construct" << nl;

    assert( new_entries.size() == counter );

    SparseMatrix ret( left.getdimout(), right.getdimin(), new_entries );
        
    // beacon.ping("Re-sort"); // LOG << "--- Sort and compress again" << nl;
    
    ret.sortandcompressentries();
    
    return ret;
    
}

// {
// 
//     LOG << "--- SparseMatrix Product" << nl;
//     LOG << "--- Sort and compress" << nl;
//     
//     left.sortandcompressentries( SparseMatrix::MatrixEntrySorting::columnwise );
//     right.sortandcompressentries( SparseMatrix::MatrixEntrySorting::rowwise );
// 
//     LOG << "--- Counting" << nl;
//     
//     int counter = 0;
//     for( SparseMatrix::MatrixEntry l : left.getentries()  )
//     for( SparseMatrix::MatrixEntry r : right.getentries() )
//         if( l.column == r.row ) 
//             counter++;
// 
//     LOG << "--- Assemble" << nl;
//     
//     std::vector<SparseMatrix::MatrixEntry> new_entries;
//     new_entries.reserve( counter );
//     for( SparseMatrix::MatrixEntry l : left.getentries()  )
//     for( SparseMatrix::MatrixEntry r : right.getentries() )
//         if( l.column == r.row ) 
//             new_entries.push_back( { l.row, r.column, l.value * r.value } );
//     
//     assert( new_entries.size() == counter );
// 
//     LOG << "--- Construct" << nl;
//     SparseMatrix ret( left.getdimout(), right.getdimin(), new_entries );
//         
//     LOG << "--- Sort and compress again" << nl;
//     ret.sortandcompressentries();
//     
//     return ret;
//     
// }





DiagonalOperator InverseDiagonalPreconditioner( const SparseMatrix& mat )
{

    assert( mat.getdimin() == mat.getdimout() );

    FloatVector diag( mat.getdimin(), 0. );

    const std::vector<SparseMatrix::MatrixEntry>& entries = mat.getentries();

    for( const auto& entry : entries )
        if( entry.row == entry.column ) 
            diag.at( entry.row ) += absolute( entry.value );

//     for( int i = 0; i < diag.getdimension(); i++ )
//         assert( diag.at( i ) > 0. );
    
    for( int i = 0; i < diag.getdimension(); i++ )
        if( not is_numerically_small( diag.at( i ) ) )
            diag.at( i ) = 1. / diag.at( i );
        else
            diag.at( i ) = 0.;
    
    return DiagonalOperator( diag );

}



Float norm_sq_of_vector( const SparseMatrix& A, const FloatVector& vec )
{
    assert( A.issquare() );
    assert( A.getdimin() == vec.getdimension() );

    const auto& entries = A.getentries();

    assert( entries.size() > 0 );

    long double ret = 0.;
    
    for( const auto& entry : entries )
    {
        // LOGPRINTF("%d %d %f\n", entry.row, entry.column, entry.value );
        ret += entry.value * vec[ entry.row ] * vec[ entry.column ];
    }
    // LOG << ret << nl;

    long double ret1 = 0.;
    long double ret2 = 0.;
    
    for( const auto& entry : entries )
    {
        // LOGPRINTF("%d %d %f\n", entry.row, entry.column, entry.value );
        if( entry.row == entry.column )
            ret1 += entry.value * vec[ entry.row ] * vec[ entry.row ];
        if ( entry.row < entry.column )
            ret2 += entry.value * vec[ entry.row ] * vec[ entry.column ];
    }
    
    // return static_cast<Float>( ret1 + 2. * ret2 );
    return static_cast<Float>( ret );
}


