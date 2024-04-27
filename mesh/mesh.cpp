

#include <cmath>
#include <algorithm>
#include <limits>


#include "../basic.hpp"
#include "../utility/random.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../dense/factorization.hpp"
#include "../dense/factorization.hpp"
#include "mesh.hpp"


#ifdef NDEBUG
const int Mesh::nullindex = std::numeric_limits<int>::max(); 
#else
const int Mesh::nullindex = -17; 
#endif


Mesh::Mesh( int inner, int outer )
: innerdimension(inner), 
  outerdimension(outer),
  coordinates(outer,0)
{
  
  /* Build up the static auxiliary data */
  /*
  auxdata.resize( innerdimension+1 );
  for( int sup = 0; sup <= innerdimension; sup++ ) 
  {
  
    auxdata[sup].resize( sup+1 );
    
    for( int sub = 0; sub <= sup; sub++ )
    {
        IndexRange from( 0, sub );
        IndexRange to( 0, sup );
        std::vector<IndexMap> sigmas = generateSigmas( from, to );
        auxdata[ sup ][ sub ] = sigmas;
//         auxdata[ std::pair<int,int>(sup,sub) ] = sigmas;
    }
  
  }
  */

  Mesh::check();
  
}

Mesh::~Mesh()
{
  Mesh::check();
}


int Mesh::getinnerdimension() const
{
  return innerdimension;
}

int Mesh::getouterdimension() const
{
  return outerdimension;
}

Coordinates& Mesh::getcoordinates()
{
  return coordinates;
}

const Coordinates& Mesh::getcoordinates() const
{
  return coordinates;
}


void Mesh::check() const 
{
  
  #if defined(DO_NOT_CHECK_MESHES)
  #warning Disabled check for Simplicial Mesh Base Class
  return;
  #endif

  #ifdef NDEBUG
  return;
  #endif
    
    // check dimension 
  assert( innerdimension >= 0 );
  assert( outerdimension >= 0 );
  assert( outerdimension == coordinates.getdimension() );
  
  // the vertices and volumes must be counted 
  // assert( has_dimension_counted( innerdimension ) );
  // assert( has_dimension_counted( 0 ) );
  
  // the counting of the vertices must agree 
  // assert( count_simplices(0) == coordinates.getnumber() );
  
  /* * Data integrity
     * 
     * - for each simplex, check that the subsimplices are non-null and unique
     * - for each supersimplex, check that the supersimplices contain that simplex
     * - morphism property of inclusion 
     *   
     */
  
  // TODO: Implement 
  
  // TODO: Develop indexing rules 
  
}

// void Mesh::print( std::ostream& out ) const
// {
//   out << text();
// }






int Mesh::index_from_pair( int sup, int sub ) const
{
  assert( 0 <= sub && sub <= sup && sup <= innerdimension );
  return sub + (sup+1) * sup / 2;
}

void Mesh::index_to_pair( int index, int& sup, int& sub ) const
{
  // FIXME: improve this incredibly ineffecient computation.
  assert( 0 <= sub && sub <= sup && sup <= innerdimension );
  for( int _sup = 1; _sup <= innerdimension; _sup++ )
  for( int _sub = 0; _sub <=            sup; _sub++ )
    if( index == index_from_pair( _sup, _sub ) )
    {
      sup = _sup; sub = _sub; return;
    }
}

int Mesh::count_subsimplices( int sup, int sub ) const 
{
  assert( 0 <= sub && sub <= sup && sup <= innerdimension );
  return binomial_integer( sup + 1, sub + 1);
}




 /*
  * 
  * Accessing subsimplices 
  *
  */

bool Mesh::is_subsimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  return im.has_value_in_range( cellsub );
}

int  Mesh::get_subsimplex_index( int sup, int sub, int cellsup, int cellsub ) const
{
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  assert( im.has_value_in_range( cellsub ) );
  return im.preimageof( cellsub );
}

int Mesh::get_subsimplex( int sup, int sub, int cellsup, int localindex ) const
{
  const IndexMap im = getsubsimplices( sup, sub, cellsup );
  return im[ localindex ];  
}

int Mesh::get_opposite_subsimplex_index( int sup, int sub, int cellsup, int localindex ) const
{
    assert( 0 <= sub && sub < sup && sup <= getinnerdimension() );
    assert( 0 <= cellsup && cellsup <= count_simplices(sup) );
    assert( 0 <= cellsup && cellsup <= count_simplices(sup) );
    
    const int cellsub = get_subsimplex( sup, sub, cellsup, localindex );
    
    assert( 0 <= cellsub && cellsub <= count_simplices(sub) );
    
    const auto my_vertices = getsubsimplices( sub, 0, cellsub );
    
    for( int opposite_index = 0; opposite_index < count_subsimplices(sup,sup-sub-1); opposite_index++ )
    {
        assert( 0 <= opposite_index && opposite_index <= count_subsimplices(sup,sup-sub-1) );
    
        const int opposite_cell = get_subsimplex( sup, sup-sub-1, cellsup, opposite_index );
        
        assert( 0 <= opposite_cell && opposite_cell <= count_simplices(sup-sub-1) );
    
        const auto other_vertices = getsubsimplices( sup-sub-1, 0, opposite_cell ); // BUG HERE????
        
        bool alive = true;
        for( int i = 0; i <    my_vertices.getSourceRange().cardinality() and alive; i++ )
        for( int j = 0; j < other_vertices.getSourceRange().cardinality() and alive; j++ )
          alive = alive && ( my_vertices[i] != other_vertices[j] );

        if( alive )
          return opposite_index;
    }
    unreachable();
}





 /*
  * 
  * Accessing supersimplices 
  *
  */

bool Mesh::is_supersimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  
  if( has_subsimplices_listed( sup, sub ) ) {
    
    return is_subsimplex( sup, sub, cellsup, cellsub );
    
  } else {
    
    assert( has_supersimplices_listed( sup, sub ) );
    std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
    return std::find( parents.begin(), parents.end(), cellsup ) != parents.end();
    
  }
  
}

int Mesh::get_firstparent_of_subsimplex( int sup, int sub, int cellsub ) const
{
  assert( has_supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  return parents[0];
}

int Mesh::get_nextparent_of_subsimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  assert( has_supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  auto it = std::find( parents.begin(), parents.end(), cellsup );
  assert( it != parents.end() );
  it++;
  if( it == parents.end() )
    return nullindex;
  else
    return *it;
}

int Mesh::get_nextparent_by_localindex( int sup, int sub, int cellsup, int localindex ) const
{
  assert( has_supersimplices_listed( sup, sub ) );
  assert( has_subsimplices_listed( sup, sub ) );
  int cellsub = getsubsimplices( sup, sub, cellsup )[localindex];
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  auto it = std::find( parents.begin(), parents.end(), cellsup );
  assert( it != parents.end() );
  it++;
  if( it == parents.end() )
    return nullindex;
  else
    return *it;
}

int Mesh::get_index_of_supersimplex( int sup, int sub, int cellsup, int cellsub ) const
{
  assert( has_supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  auto it = std::find( parents.begin(), parents.end(), cellsup );
  assert( it != parents.end() );
  return it - parents.begin();
}

int Mesh::get_supersimplex_by_index( int sup, int sub, int cellsub, int parentindex ) const
{
  assert( has_supersimplices_listed( sup, sub ) );
  std::vector<int> parents = getsupersimplices( sup, sub, cellsub );
  return parents[ parentindex ];
}













void Mesh::set_flags( int dim, SimplexFlag flag )
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( has_dimension_counted( dim ) );
    for( int s = 0; s < count_simplices(dim); s++ )
        set_flag( dim, s, flag );
}

const std::vector<SimplexFlag> Mesh::get_flags( int dim ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( has_dimension_counted( dim ) );
    std::vector<SimplexFlag> flags( count_simplices(dim), SimplexFlag::SimplexFlagInvalid );
    for( int s = 0; s < count_simplices(dim); s++ )
        flags[s] = get_flag( dim, s );
    return flags;
}

void Mesh::set_flags( int dim, std::vector<SimplexFlag> flags )
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( has_dimension_counted( dim ) );
    assert( flags.size() == count_simplices( dim ) );
    for( int s = 0; s < count_simplices(dim); s++ )
        set_flag( dim, s, flags[s] );
}

void Mesh::automatic_dirichlet_flags()
{
    const int full = getinnerdimension();
    
    assert( has_dimension_counted(full-1) );
    assert( has_supersimplices_listed(full,full-1) );
    
    for( int d = 0; d <= getinnerdimension(); d++ )
        set_flags( d, SimplexFlag::SimplexFlagNull );
    
    for( int s = 0; s < count_simplices(full-1); s++ )
        if( 
            get_firstparent_of_subsimplex( full, full-1, s ) == nullindex 
            || 
            get_nextparent_of_subsimplex( full, full-1, get_firstparent_of_subsimplex( full, full-1, s ), s ) == nullindex
        ) {
            set_flag( full-1, s, SimplexFlag::SimplexFlagDirichlet );
        }
        
    for( int s = 0; s < count_simplices(full-1); s++ )
        if( get_flag( full-1, s ) == SimplexFlag::SimplexFlagDirichlet )
            for( int d = 0; d < full-1; d++ )
                for( int subindex = 0; subindex < count_subsimplices( full-1, d ); subindex++ ) 
                    set_flag( d, get_subsimplex( full-1, d, s, subindex ), SimplexFlag::SimplexFlagDirichlet );

}


void Mesh::check_dirichlet_flags( bool check_for_full_dirichlet )
{
    const int full = getinnerdimension();
    
    assert( has_dimension_counted(full-1) );
    assert( has_supersimplices_listed(full,full-1) );
    
    // check that full-dimensional simplices have no Dirichlet condition 
    for( int s = 0; s < count_simplices(full); s++ )
        assert( get_flag( full, s ) == SimplexFlag::SimplexFlagNull );

    // check that no flags are invalid
    for( int d = 0; d <= full-1; d++ )
    for( int s = 0; s < count_simplices(d); s++ )
        assert( get_flag( d, s ) != SimplexFlag::SimplexFlagInvalid );
    
    // if some cell has the Dirichlet condition, then so do its subcells
    for( int d = 1; d <= full-1; d++ )
    for( int s = 0; s < count_simplices(d); s++ )
        if( get_flag( d, s ) == SimplexFlag::SimplexFlagDirichlet )
            for( int subindex = 0; subindex < count_subsimplices( d, d-1 ); subindex++ ) 
                assert( get_flag( d-1, get_subsimplex( d, d-1, s, subindex ) ) == SimplexFlag::SimplexFlagDirichlet );
    
    // for all lower dimensional simplices 
    // if they are Dirichlet, check that they are contained in a Dirichlet face 
    // if they are not, then check that they are not contained in a Dirichlet face 
    for( int d = 0; d <= full-2; d++ )
    for( int sub = 0; sub < count_simplices(d); sub++ ) 
    {

        if( get_flag( d, sub ) == SimplexFlag::SimplexFlagDirichlet ) {

            bool found = false;
            
            for( 
                int sup = get_firstparent_of_subsimplex( d+1, d, sub ); 
                sup != nullindex; 
                sup = get_nextparent_of_subsimplex( d+1, d, sup, sub ) 
            ){
                found = found or ( get_flag(d+1,sup) == SimplexFlag::SimplexFlagDirichlet );
            }
            
            if( not found ) LOG << d << space << sub;
            
            assert( found );

        } else {

            assert( get_flag( d, sub ) == SimplexFlag::SimplexFlagNull );

            for( 
                int sup = get_firstparent_of_subsimplex( d+1, d, sub ); 
                sup != nullindex; 
                sup = get_nextparent_of_subsimplex( d+1, d, sup, sub ) 
            ){
                assert( get_flag(d+1,sup) == SimplexFlag::SimplexFlagNull ); // turn into equality check
            }
        }

    }

        
    if( not check_for_full_dirichlet ) return;


    
    // for each face, check that those faces with less than two parents have the Dirichlet flag 
    for( int s = 0; s < count_simplices(full-1); s++ )
        if( get_firstparent_of_subsimplex( full, full-1, s ) == nullindex 
            || 
            get_nextparent_of_subsimplex( full, full-1, get_firstparent_of_subsimplex( full, full-1, s ), s ) == nullindex ) {
            assert( get_flag( full-1, s ) == SimplexFlag::SimplexFlagDirichlet );
        } else {
            assert( get_flag( full-1, s ) == SimplexFlag::SimplexFlagNull      );
        }
 
    // for each face, check that those faces with at least two parents have no Dirichlet flag 
    for( int s = 0; s < count_simplices(full-1); s++ )
        if( get_firstparent_of_subsimplex( full, full-1, s ) != nullindex 
            && 
            get_nextparent_of_subsimplex( full, full-1, get_firstparent_of_subsimplex( full, full-1, s ), s ) != nullindex
        ) {
            assert( get_flag( full-1, s ) == SimplexFlag::SimplexFlagNull      );
        } else {
            assert( get_flag( full-1, s ) == SimplexFlag::SimplexFlagDirichlet );
        }
 
    
    
}















Float Mesh::getDiameter( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix vcm = getVertexCoordinateMatrix( dim, index );
    
    DenseMatrix dist(dim+1,dim+1);
    for( int i = 0; i <= dim; i++ )
    for( int j = 0; j <= dim; j++ )
        dist(i,j) = ( vcm.getcolumn(i) - vcm.getcolumn(j) ).norm();
    
    assert( dist.isnonnegative() && dist.isfinite() );
    
    return dist.maxabsoluteentry();
}

Float Mesh::getMaximumDiameter() const
{
    Float ret = 0.;
    for( int e = 0; e < count_simplices(1); e++ )
        ret = maximum( ret, getDiameter(1,e) );
    return ret;
}

Float Mesh::getMinimumDiameter() const
{
    Float ret = std::numeric_limits<Float>::infinity();
    for( int e = 0; e < count_simplices(1); e++ )
        ret = minimum( ret, getDiameter(1,e) );
    assert( std::isfinite(ret) );
    return ret;
}
        
        

Float Mesh::getMeasure( int dim, int index ) const 
{
    assert( 0 <= dim   && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix Jac = getTransformationJacobian( dim, index );
    
//     DenseMatrix temp = Transpose( Jac ) * Jac;

    Float det;
    
    if( Jac.issquare() )
        det = absolute( Determinant( Jac ) ); 
    else
        det = std::sqrt( absolute( Determinant( Transpose(Jac) * Jac ) ) );

    return det / factorial_numerical( dim );
}

Float Mesh::getHeight( int dim, int cell, int vertexindex ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= cell && cell <= count_simplices(dim) );
    assert( 0 <= vertexindex && vertexindex <= dim );
    
    int oppositeface_index = get_opposite_subsimplex_index( dim, 0, cell, vertexindex );
    int oppositeface = get_subsimplex( dim, dim-1, cell, oppositeface_index );
    
    Float vol_t = getMeasure( dim, cell );
    Float vol_f = getMeasure( dim-1, oppositeface );

    // vol_t = vol_f * height / dim 

    return ( vol_t / vol_f ) * (dim);
}


Float Mesh::getHeightQuotient( int dim, int cell ) const
{
    assert( 1 <= dim && dim <= getinnerdimension() );
    assert( 0 <= cell && cell <= count_simplices(dim) );
    
    Float vol_t = getMeasure( dim, cell );
    Float diameter = getDiameter( dim, cell );
    Float ret = 0.;
    for( int i = 0; i <= dim; i++ )
    {
        Float vol_f = getMeasure( dim-1, get_subsimplex(dim,dim-1,cell,i) );
        Float height = ( vol_t / vol_f ) * (dim);
        Float ratio = diameter / height;
        ret = maximum( ret, ratio );
    }
    
    return ret;
}

Float Mesh::getHeightQuotient( int dim ) const
{
    assert( 1 <= dim && dim <= getinnerdimension() ); 
    
    Float height_ratio = 0.;
    for( int s = 0; s < count_simplices(dim); s++ )
        height_ratio = maximum( height_ratio, getHeightQuotient(dim,s) );
    
    return height_ratio;
}

Float Mesh::getHeightQuotient() const
{
    Float ret = 0.;
    for( int d = 1; d <= getinnerdimension(); d++ )
        ret = maximum( ret, getHeightQuotient(d) ); 
    return ret;
}
        
        

Float Mesh::getShapemeasure( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    return power_numerical( getDiameter( dim, index ), dim ) / getMeasure( dim, index );
}

Float Mesh::getShapemeasure( int dim ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    
    Float shapemeas = 0.;
    for( int s = 0; s < count_simplices(dim); s++ )
        shapemeas = maximum( shapemeas, power_numerical( getDiameter( dim, s ), dim ) / getMeasure( dim, s ) );
    return shapemeas;
}

Float Mesh::getShapemeasure() const 
{
    return getShapemeasure( getinnerdimension() ); 
}


int Mesh::getVertexPatchSize() const
{
    int ret = 0;
    for( int s = 0; s < count_simplices(0); s++ )
    {
        int count = getsupersimplices( getinnerdimension(), 0, s ).size();
        ret = maximum(ret,count);
    }
    return ret;
}

int Mesh::getSupersimplexSize( int dim ) const
{
    assert( 0 <= dim && dim < getinnerdimension() );
    
    int ret = 0;
    for( int s = 0; s < count_simplices(dim); s++ )
    {
        int count = getsupersimplices( dim+1, dim, s ).size();
        ret = maximum(ret,count);
    }
    return ret;
}
        
        

Float Mesh::getComparisonQuotient() const
{
    int n = getinnerdimension();
    
    Float physical_ratio = 1.;
    for( int s = 0; s < count_simplices(n); s++ )
    {
        auto edges = getsubsimplices( n, 1, s ).getvalues();
        Float diameter = getDiameter( n, s );
        for( auto e : edges )
        {
            Float edgelength = getDiameter( 1, e );
            physical_ratio = maximum( physical_ratio, diameter / edgelength );
        }
    }

    Float vertex_ratio = 1.;
    for( int v = 0; v < count_simplices(0); v++ )
    {
        auto volumes = getsupersimplices( n, 0, v );
        auto edges  = getsupersimplices( 1, 0, v );

        Float max_diameter = 0.;
        for( auto t : volumes ) max_diameter = maximum( max_diameter, getDiameter( n, t ) );
        
        Float min_length = 0.;
        for( auto e : edges ) min_length = minimum( min_length, getDiameter( 1, e ) );

        assert( max_diameter > min_length );

    }
    
    return maximum( physical_ratio, vertex_ratio );

}

Float Mesh::getRadiiQuotient( int dim ) const
{
    int n = getinnerdimension();

    return getHeightQuotient() * (dim+1);
    
    Float comparison_quotient = getComparisonQuotient();
    
    Float shape_measure = getShapemeasure();
    
    Float sigma = factorial_numerical(n) / ( power_numerical( n, n / 2. ) * shape_measure * std::sqrt(2) );
    
    return comparison_quotient * std::sqrt(n) * (n+1) / sigma;
}





FloatVector Mesh::get_midpoint( int dim, int index ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    FloatVector mid( getouterdimension(), 0. );
    
    for( int v = 0; v <= dim; v++ )
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] += getcoordinates().getdata( get_subsimplex( dim, 0, index, v ), d );
    
    for( int d = 0; d < getouterdimension(); d++ )
      mid[d] /= dim + 1;

    return mid;
}

FloatVector Mesh::getPointFromBarycentric( int dim, int index, const FloatVector& barycoords ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    assert( barycoords.getdimension() == dim+1 );
    assert( is_numerically_close( barycoords.sum(), 1. ) );

    FloatVector ret( getouterdimension(), 0. );
    
    for( int v = 0; v <= dim; v++ )
    for( int d = 0; d < getouterdimension(); d++ )
      ret[d] += barycoords[v] * getcoordinates().getdata( get_subsimplex( dim, 0, index, v ), d );
    
    return ret;
}

FloatVector Mesh::get_random_point( int dim, int index ) const
{
    auto randomcoords = get_random_barycentric_coordinates(dim);
    return getPointFromBarycentric( dim, index, randomcoords );
    
    // std::vector<Float> barycoords( dim+2, 0. );
    // for( auto& x : barycoords ) x = random_uniform();
    // barycoords.back() = 1.;
    // barycoords.front() = 0.;
    // for( auto& x : barycoords ) assert( 0. <= x and x <= 1. );
    // std::sort( barycoords.begin(), barycoords.end() );
    // for( int i = 0; i <= dim; i++ ) barycoords[i] = barycoords[i+1] - barycoords[i];
    // barycoords.pop_back();
    // return getPointFromBarycentric( dim, index, FloatVector(barycoords) );
}




int Mesh::get_longest_edge_index( int dim, int index ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    const int N = binomial_integer(dim+1,2);
    
    int longest_edge_index = 0;

    Float diameter = getDiameter(1,0);

    for( int n = 1; n < N; n++ )
    {
        int e = get_subsimplex( dim, 1, index, n );
        Float diameter_e = getDiameter(1,e);
        if( diameter_e > diameter ){
            diameter = diameter_e;
            longest_edge_index = e;
        }
    }
    
    return longest_edge_index;
}
        
        

DenseMatrix Mesh::getVertexCoordinateMatrix( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix ret( getouterdimension(), dim+1 );
    
    for( int v = 0; v <= dim; v++ )
    for( int d = 0; d < getouterdimension(); d++ )
        ret( d, v ) = coordinates.getdata( get_subsimplex( dim, 0, index, v ), d );
    
    return ret;
}


DenseMatrix Mesh::getTransformationJacobian( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix ret( getouterdimension(), dim );
    
    DenseMatrix vcm = getVertexCoordinateMatrix( dim, index );
    
    for( int v = 0; v < dim; v++ )
    for( int c = 0; c < getouterdimension(); c++ )
        ret( c, v ) = vcm( c, v+1 ) - vcm( c, 0 );
    
    return ret;
}


DenseMatrix Mesh::getBarycentricProjectionMatrix( int dim, int index ) const
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    const auto Jac = getTransformationJacobian( dim, index );

    assert( Jac.getdimout() >= Jac.getdimin() );
    
    DenseMatrix ret( Jac.getdimin()+1, Jac.getdimout(), 0.0 );
    
    for( int r = 0; r < Jac.getdimin();  r++ )
    for( int c = 0; c < Jac.getdimout(); c++ )
        ret( r+1, c ) = Jac(c,r);
    
    assert( ret.isfinite() );
    
    return ret;
}


// DenseMatrix Mesh::getGradientMatrix( int dim, int index ) const 
// {
//     assert( 0 <= dim && dim <= getinnerdimension() );
//     assert( 0 <= index && index < count_simplices(dim) );
    
//     DenseMatrix ret( getouterdimension(), dim+1, 0. );
    
//     DenseMatrix vcm = getVertexCoordinateMatrix( dim, index );
    
//     for( int v = 1; v <= dim; v++ )
//     for( int c = 0; c < getouterdimension(); c++ )
//     {
//         ret( c, v ) = vcm( c, v ) - vcm( c, 0 );
//         ret( c, 0 ) -= ret( c, v );
//     }
    
//     return ret;
// }
DenseMatrix Mesh::getGradientMatrix( int dim, int index ) const 
{
    assert( 0 <= dim   && dim   <= getinnerdimension() );
    assert( 0 <= index && index <  count_simplices(dim) );
    
    DenseMatrix multiplier( dim, dim+1, 0. );
    for( int i = 0; i <  dim; i++ ) {
        multiplier(i,i+1) =  1.;
        multiplier(i,  0) = -1.;
    }
    
    // D^-1 D^-t = ( D^t D )^-1
    DenseMatrix Jac    = getTransformationJacobian( dim, index );

    DenseMatrix R( Jac.getdimin() );
    DenseMatrix Q( Jac.getdimout(), Jac.getdimin() );
    QRFactorization( Jac, Q, R );

    return Q * Transpose( Inverse(R) ) * multiplier;
}
        
DenseMatrix Mesh::getGradientProductMatrix( int dim, int index ) const 
{
    assert( 0 <= dim   && dim   <= getinnerdimension() );
    assert( 0 <= index && index <  count_simplices(dim) );
    
    DenseMatrix multiplier( dim, dim+1, 0. );
    for( int i = 0; i <  dim; i++ ) {
        multiplier(i,i+1) =  1.;
        multiplier(i,  0) = -1.;
    }
    
    // D^-1 D^-t = ( D^t D )^-1
    
    DenseMatrix Jac    = getTransformationJacobian( dim, index );
    
    // DenseMatrix middle = Inverse( Transpose(Jac) * Jac );
    // return Transpose(multiplier) * middle * multiplier;

    DenseMatrix R( Jac.getdimin() );
    DenseMatrix Q( Jac.getdimout(), Jac.getdimin() );
    QRFactorization( Jac, Q, R );
    DenseMatrix Rinv( Inverse(R) );
    return Transpose(multiplier) * ( Rinv * Transpose(Rinv) ) * multiplier;
}
        
DenseMatrix Mesh::getGradientProductMatrixRightFactor( int dim, int index ) const 
{
    assert( 0 <= dim && dim <= getinnerdimension() );
    assert( 0 <= index && index < count_simplices(dim) );
    
    DenseMatrix multiplier( dim, dim+1, 0. );
    for( int i = 0; i <  dim; i++ ) {
        multiplier(i,i+1) =  1.;
        multiplier(i,  0) = -1.;
    }
    
    // D^-1 D^-t = ( D^t D )^-1
    DenseMatrix Jac    = getTransformationJacobian( dim, index );
    DenseMatrix middle = Inverse( Transpose(Jac) * Jac ); //Transpose(Jac) * Jac; //TODO: Understand
    
    DenseMatrix middle_rightfactor = Transpose( CholeskyDecomposition( middle ) ); 

    assert( ( middle - Transpose(middle_rightfactor) * middle_rightfactor ).is_numerically_small() );
    
    return middle_rightfactor * multiplier; //TODO Probelesen
}


FloatVector get_random_barycentric_coordinates( int dim )
{
    FloatVector samples( dim );
    do samples.random_within_range(0.,1.); while( samples.sum() > 1. );
    FloatVector randomcoords( dim+1 );
    for( int p = 0; p < dim; p++ ) randomcoords[p] = samples[p];
    randomcoords[dim] = 1. - samples.sum();
    return randomcoords;
}   











FloatVector Mesh::transform_whitney_to_euclidean( int dim, const FloatVector& whitneyvalues, int zero_padding ) const 
{
    Assert( whitneyvalues.getdimension() == count_simplices(dim) * (dim+1) );

    FloatVector ret( ( getouterdimension() + zero_padding ) * count_simplices(dim), 0. );

    for( int t = 0; t < count_simplices(dim); t++ )
    {
        FloatVector coefficients( dim+1 );
        
        for( int i = 0; i <= getinnerdimension(); i++ )
            coefficients[i] = whitneyvalues.at( t * (dim+1) + i );
        
        FloatVector directions = getGradientMatrix(dim,t) * coefficients;

        ret.setslice( ( getouterdimension() + zero_padding ) * t, directions );
    }
        
    return ret;
}






void Mesh::shake_interior_vertices( Float intensity, Float probability )
{
    assert( 0. <= probability and probability <= 1.0 );
    assert( 0. <= intensity and intensity <= 1.0 );
    assert( getinnerdimension() == getouterdimension() ); 
    
    const int dim = getinnerdimension();

    for( int v = 0; v < count_simplices(0); v++ )
    {
        if( random_uniform() > probability ) continue; 

        const auto es = getsupersimplices( 1, 0, v );

        bool boundary_detected = false;

        for( const int e : es )
        {
            const auto ts = getsupersimplices( dim, 1, e );
            assert( ts.size() >= 1 );
            if( ts.size() == 1 ) boundary_detected = true;
        }

        if( boundary_detected ) 
            continue;

        Float radius = std::numeric_limits<Float>::infinity();

        const auto ts = getsupersimplices( dim, 0, v );

        for( const int t : ts )
        {
            int vi = get_subsimplex_index( dim, 0, t, v );
            Float height = getHeight( dim, t, vi );
            radius = minimum( height, radius );
        }

        assert( std::isfinite(radius) && radius >= 0. );

        FloatVector shift(dim,0.);
        shift.random_within_range(-1.,1.);
        shift /= shift.l2norm();
        assert( shift.isfinite() );
        shift *= sqrt( random_uniform() ) * radius * intensity;

        for( int c = 0; c < dim; c++ )
            getcoordinates().setdata( v, c, getcoordinates().getdata( v, c ) + shift[c] );

    }
}
