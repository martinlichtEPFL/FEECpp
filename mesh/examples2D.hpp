#ifndef INCLUDEGUARD_EXAMPLES_2D_HPP
#define INCLUDEGUARD_EXAMPLES_2D_HPP


#include <cmath>
#include <vector>

#include "mesh.hpp"
#include "mesh.simplicial2D.hpp"











inline MeshSimplicial2D StandardSquare2D_simple()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 4, {
        -1., -1., // 0
        -1.,  1., // 1
         1., -1., // 2
         1.,  1.  // 3
      } ),
      {
        { 0, 1, 3 },
        { 0, 2, 3 }
      }
    );
}


inline MeshSimplicial2D StandardSquare2D_alternative()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 4, {
        -1.,  1., // 0
         1., -1., // 1
        -1., -1., // 2
         1.,  1.  // 3
      } ),
      {
        { 0, 1, 2 },
        { 0, 1, 3 }
      }
    );
}


inline MeshSimplicial2D StandardSquare2D_tiles3x3()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 16, {
         0.,  0., // 0
         0.,  1., // 1
         0.,  2., // 2
         0.,  3., // 3
         1.,  0., // 4
         1.,  1., // 5
         1.,  2., // 6
         1.,  3., // 7
         2.,  0., // 8
         2.,  1., // 9
         2.,  2., // 10
         2.,  3., // 11
         3.,  0., // 12
         3.,  1., // 13
         3.,  2., // 14
         3.,  3., // 15
      } ),
      {
        // Ecken
        { 0, 1, 4 },
        { 1, 4, 5 },

        { 2, 3, 7 },
        { 2, 6, 7 },
        
        { 8, 9,12 },
        { 9,12,13 },
        
        {10,11,14 },
        {11,14,15 },
        // Horizontal, vertikal
        { 1, 2, 5 },
        { 2, 5, 6 },
        
        { 6, 7,10 },
        { 7,10,11 },

        { 4, 5, 8 },
        { 5, 8, 9 },

        { 9,10,13 },
        {10,13,14 },
        // mitte 
        { 5, 6, 9 },
        { 6, 9,10 }
      }
    );
}


inline MeshSimplicial2D StandardSquare2D_strange14()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 12, {
        -1.0, -1.0, // 0
        -1.0,  1.0, // 1
         1.0, -1.0, // 2
         1.0,  1.0, // 3
         //
        -1.0, 0.1, // 4
        -0.2,-1.0, // 5
         0.3, 1.0,  // 6
         1.0,-0.2,  // 7
         //
        -0.3, -0.5, // 8
        -0.4,  0.4, // 9
         0.4, -0.5, // A
         0.5,  0.3, // B
         //
        
      } ),
      {
        {  0,  4,  8 },
        {  0,  5,  8 },
        {  5,  8, 10 },
        {  2,  5, 10 },
        {  2,  7, 10 },
        
        {  4,  8,  9 },
        {  8,  9, 11 },
        {  8, 10, 11 },
        {  7, 10, 11 },
        
        {  1,  4,  9 },
        {  1,  6,  9 },
        {  6,  9, 11 },
        {  3,  6, 11 },
        {  3,  7, 11 } 
      }
    );
}


inline MeshSimplicial2D StandardSquare2D_centered()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 9, {
         0.,  0., // 0
         1.,  0., // 1
         1.,  1., // 2
         0.,  1., // 3
        -1.,  1., // 4
        -1.,  0., // 5
        -1., -1., // 6
         0., -1., // 7
         1., -1.  // 8
      } ),
      {
        { 0, 1, 2 },
        { 0, 2, 3 },
        { 0, 3, 4 },
        { 0, 4, 5 },
        { 0, 5, 6 },
        { 0, 6, 7 },
        { 0, 7, 8 },
        { 0, 8, 1 }
      }
    );
}


inline MeshSimplicial2D StandardSquare2D()
{
    return StandardSquare2D_strange14();
}






inline MeshSimplicial2D UnitSquare2D()
{
    auto M = StandardSquare2D();
    
    M.getcoordinates().shift( FloatVector{1.,1.} );
    M.getcoordinates().scale( 0.5 );
    
    return M;
}

inline MeshSimplicial2D UnitSquare2D_alternative()
{
    auto M = StandardSquare2D_alternative();
    
    M.getcoordinates().shift( FloatVector{1.,1.} );
    M.getcoordinates().scale( 0.5 );
    
    return M;
}

inline MeshSimplicial2D UnitSquare2D_centered()
{
    auto M = StandardSquare2D_centered();
    
    M.getcoordinates().shift( FloatVector{1.,1.} );
    M.getcoordinates().scale( 0.5 );
    
    return M;
}

inline MeshSimplicial2D UnitSquare2D_strange14()
{
    auto M = StandardSquare2D_strange14();
    
    M.getcoordinates().shift( FloatVector{1.,1.} );
    M.getcoordinates().scale( 0.5 );
    
    return M;
}
















inline MeshSimplicial2D Triangle_by_angle_and_length( Float angle, Float length )
{
    assert( 0 <= angle && angle <= Constants::twopi );
    assert( 0 <= length );
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 3, {
        0., 0., // 0
        1., 0., // 1
        length * cos(angle) , length * sin(angle)  // 2
      } ),
      {
        { 0, 1, 2 }
      }
    );
}



inline MeshSimplicial2D UnitTriangle2D()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 3, {
        0., 0., // 0
        1., 0., // 1
        0., 1.  // 2
      } ),
      {
        { 0, 1, 2 }
      }
    );
}





inline MeshSimplicial2D RegularSimplex2D()
{
    return MeshSimplicial2D(
      3,  
      Coordinates( 3, 3, {
         0., 0., 0., // 0
         1., 1., 0., // 1
         1., 0., 1., // 2
      } ),
      {
        { 0, 1, 2 }
      }
    );
}




inline MeshSimplicial2D TetrahedralSurface2D()
{
    return MeshSimplicial2D(
      3,
      Coordinates( 3, 4, {
         0.,  0.,  0., // 0
         1.,  0.,  0., // 1
         0.,  1.,  0., // 2
         0.,  0.,  1.,  // 3
      } ),
      {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 }
      }
    );
}




// inline MeshSimplicial2D LShapedDomain2D()
// {
//     return MeshSimplicial2D(
//       2,
//       Coordinates( 2, 8, {
//          0.,  0.,  // 0
//          1.,  0.,  // 1
//          1.,  1.,  // 2
//          0.,  1.,  // 3
//         -1.,  1.,  // 4
//         -1.,  0.,  // 5
//         -1., -1.,  // 6
//          0., -1.   // 7
//       } ),
//       {
//         { 0, 1, 2 },
//         { 0, 2, 3 },
//         { 0, 3, 4 },
//         { 0, 4, 5 },
//         { 0, 5, 6 },
//         { 0, 6, 7 }
//       }
//     );
// }
// inline MeshSimplicial2D LShapedDomain2D()
// {
//     return MeshSimplicial2D(
//       2,
//       Coordinates( 2, 8, {
//          0.,  0.,  // 0
//          0., -1.,  // 1
//         -1., -1.,  // 2
//         -1.,  0.,  // 3
//         -1.,  1.,  // 4
//          0.,  1.,  // 5
//          1.,  1.,  // 6
//          1.,  0.   // 7
//       } ),
//       {
//         { 0, 1, 2 },
//         { 0, 2, 3 },
//         { 0, 3, 4 },
//         { 0, 4, 5 },
//         { 0, 5, 6 },
//         { 0, 6, 7 }
//       }
//     );
// }
inline MeshSimplicial2D LShapedDomain2D()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 6, {
         0.,  0.,  // 0
         0., -1.,  // 1
        -1., -1.,  // 2
        -1.,  1.,  // 3
         1.,  1.,  // 4
         1.,  0.   // 5
      } ),
      {
        { 0, 3, 2 },
        { 0, 3, 4 },
        { 0, 2, 1 },
        { 0, 4, 5 }         
      }
    );
}





inline MeshSimplicial2D SlitDomain2D()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 10, {
         0.,  0., // 0
         1.,  0., // 1 +++
         1.,  1., // 2
         0.,  1., // 3
        -1.,  1., // 4
        -1.,  0., // 5
        -1., -1., // 6
         0., -1., // 7
         1., -1., // 8
         1.,  0.  // 9 ---
      } ),
      {
        { 0, 1, 2 }, // 0
        { 0, 2, 3 }, // 1
        { 0, 3, 4 }, // 2
        { 0, 4, 5 }, // 3
        { 0, 5, 6 }, // 4
        { 0, 6, 7 }, // 5
        { 0, 7, 8 }, // 6
        { 0, 8, 9 }  // 7
      }
    );
}






inline MeshSimplicial2D RhombicAnnulus2D()
{
    return MeshSimplicial2D(
      2,
      Coordinates( 2, 8, {
         0.5,  0.0,  // 0
         1.0,  0.0,  // 1
         0.0,  0.5,  // 2
         0.0,  1.0,  // 3
        -0.5,  0.0,  // 4
        -1.0,  0.0,  // 5
         0.0, -0.5,  // 6
         0.0, -1.0   // 7
      } ),
      {
        { 0, 1, 2 },
        { 1, 2, 3 },
        { 2, 3, 4 },
        { 3, 4, 5 },
        { 4, 5, 6 },
        { 5, 6, 7 },
        { 0, 6, 7 },
        { 0, 1, 7 }
      }
    );
}








inline MeshSimplicial2D Halo( int Levels = 1, int Division = 3 )
{
    assert( Levels >= 1 );
    assert( Division >= 3 );
    
    int num_vertices = Division * (Levels+1);
    
    int num_triangles = 2 * Division * Levels;
    
    
    
    // 2. Create the coordinates 
    
    std::vector<Float> coords;
    coords.reserve( 3 * num_vertices );
    
    for( int level    = 0; level    <= Levels;   level++    )
    for( int division = 0; division <  Division; division++ ) {
            coords.push_back( std::sin( 2 * Constants::pi * division / (Float)Division ) );
            coords.push_back( std::cos( 2 * Constants::pi * division / (Float)Division ) );
            coords.push_back( 0.0 + level / (Float) Levels );
    }
    
    assert( coords.size() == 3 * num_vertices );
    
    std::vector<std::array<int,3>> tris;
    tris.reserve( num_triangles );
    
    for( int level = 0; level < Levels; level++ ) 
    for( int division = 0; division < Division; division++ )
    {
        
        tris.push_back( { 
            (level+0)*Division + (division    ) % Division,
            (level+1)*Division + (division    ) % Division, 
            (level+1)*Division + (division + 1) % Division
            } );
        
        tris.push_back( { 
            (level+0)*Division + (division    ) % Division, 
            (level+0)*Division + (division + 1) % Division,
            (level+1)*Division + (division + 1) % Division
            } );
                
    }
    
    assert( tris.size() == num_triangles );

    for( auto& t : tris ) 
    {
      if( t[0] > t[1] ) { int i = t[0]; t[0] = t[1]; t[1] = i; }
      if( t[0] > t[2] ) { int i = t[0]; t[0] = t[2]; t[2] = i; }
      if( t[1] > t[2] ) { int i = t[1]; t[1] = t[2]; t[2] = i; }
      assert( t[0] < t[1] and t[1] < t[2] );
    }
    
    return MeshSimplicial2D(
      3,
      Coordinates( 3, num_vertices, coords ),
      tris
    );
}









inline MeshSimplicial2D UnitDisk( int L = 1 )
{
    assert( L >= 1 );
    
    // 1. Calculate the number of vertices and triangles 
    
    //  3 + 6 + 12 + 24 + ... 
    //= 3 ( 1 + 2 + 4 + 8 + .... )
    //= 3 ( 2^( L+1 ) - 1 )
    int num_vertices = 3 * ((1<<L)-1);
    
    //  1 + (3+6) + (6+12) + (12+24) + ...
    //= 1 +   9   +    18  +     36  + ...
    int num_triangles = 9 * (1<<(L-1)) - 8;
    
    
    
    // 2. Create the coordinates 
    
    std::vector<Float> coords;
    coords.reserve( 2 * num_vertices );
    
    // 2.1 Calculate the radii
    
    std::vector<Float> Rs(L);
    
    Rs[0] = 1.;
    for( int l = 1; l < L; l++ ) Rs[l] = Rs[l-1] + 1.5 * Constants::twopi * Rs[l-1] / ( 3. * power_numerical( 2., l ) );
    
    Float Rmax = Rs[0];
    for( const Float R : Rs ) if( Rmax < R ) Rmax = R; // *std::max_element(Rs.begin(),Rs.end());
    
    for( Float& R : Rs ) R /= Rmax;
    
    // 2.2 fill in the values
    
    for( int l = 1; l <= L; l++ ) {
        
        int N = 3 * (1<<(l-1));
        
        Float radius = Rs[l-1];
        
        for( int a = 0; a < N; a++ )
        {
            coords.push_back( radius * std::cos( 2 * Constants::pi * a / (Float)N ) );
            coords.push_back( radius * std::sin( 2 * Constants::pi * a / (Float)N ) );
        }
    }
    
    LOG << coords.size() / 2 << space << num_vertices << nl;
    assert( coords.size() == 2 * num_vertices );
    
    std::vector<std::array<int,3>> tris;
    tris.reserve( num_triangles );
    
    tris.push_back( {0,1,2} );
    
    for( int l = 1; l < L; l++ )
    {
        int base_inner  = 3 * ( power_of_two( l-1 ) - 1 );
        int count_inner = 3 * power_of_two( l-1 );
        
        int base_outer  = base_inner + count_inner;
        int count_outer = 2 * count_inner;
        
        for( int i = 0; i < count_inner; i++ )
        {
            
            tris.push_back( { 
                base_inner + i, 
                base_outer + 2*i,
                base_outer + (2*i + 1)%count_outer,
                } );
            tris.push_back( { 
                base_inner + i, 
                base_inner + (i+1)%count_inner,
                base_outer + 2*i + 1,
                } );
            tris.push_back( { 
                base_outer + (2*i+2)%count_outer,
                base_inner + (i+1)%count_inner,
                base_outer + 2*i + 1,
                } );
            
        }
    }
    
    LOG << coords.size() / 2 << space << num_vertices << nl;
    LOG << tris.size() << space << num_triangles << nl;
    assert( tris.size() == num_triangles );
    
//     for( auto t : tris ){ for( auto v : t ) 
//         LOG << v << space; LOG << nl; }
    
    for( auto& t : tris ) //std::sort( t.begin(), t.end() );
    {
      if( t[0] > t[1] ) { int i = t[0]; t[0] = t[1]; t[1] = i; }
      if( t[0] > t[2] ) { int i = t[0]; t[0] = t[2]; t[2] = i; }
      if( t[1] > t[2] ) { int i = t[1]; t[1] = t[2]; t[2] = i; }
      assert( t[0] < t[1] and t[1] < t[2] );
    }
    
//     for( auto t : tris ){ for( auto v : t ) 
//         LOG << v << space; LOG << nl; }
    
    return MeshSimplicial2D(
      2,
      Coordinates( 2, num_vertices, coords ),
      tris
    );
}









inline MeshSimplicial2D Annulus( int Linner, int Louter = 1 )
{
    assert( Linner >= 1 && Louter > Linner );
    
    // 1. Calculate the number of vertices and triangles 
    
    //  3 + 6 + 12 + 24 + ... 
    //= 3 ( 1 + 2 + 4 + 8 + .... )
    //= 3 ( 2^( L+1 ) - 1 )
    int num_vertices  = 3 * ((1<<Louter)-1)     - ( 3 * ((1<<(Linner-1))-1) );
    
    //  1 + (3+6) + (6+12) + (12+24) + ...
    //= 1 +   9   +    18  +     36  + ...
    int num_triangles = 9 * (1<<(Louter-1)) - 8 - ( 9 * (1<<(Linner-1)) - 8 );
    
    
    
    // 2. Create the coordinates 
    
    std::vector<Float> coords;
    coords.reserve( 2 * num_vertices );
    
    // 2.1 Calculate the radii
    
    std::vector<Float> Rs(Louter);
    for( int l = 0; l < Linner; l++ ) Rs[l] = 1.;
    for( int l = Linner; l < Louter; l++ ) Rs[l] = Rs[l-1] + 1.5 * 2 * Constants::pi * Rs[l-1] / ( 3. * (1<<(l)) );
    
    Float Rmax = Rs[0];
    for( const Float R : Rs ) if( Rmax < R ) Rmax = R; // *std::max_element(Rs.begin(),Rs.end());
    
    for( Float& R : Rs ) R = power_numerical( R / Rmax, 1.5 );
    
    // 2.2 fill in the values
    
    for( int l = Linner; l <= Louter; l++ ) {
        
        int N = 3 * (1<<(l-1));
        
        Float radius = Rs[l-1];
        
        for( int a = 0; a < N; a++ )
        {
            coords.push_back( radius * std::cos( 2 * Constants::pi * a / (Float)N ) );
            coords.push_back( radius * std::sin( 2 * Constants::pi * a / (Float)N ) );
        }
    }
    
    LOG << coords.size() / 2 << space << num_vertices << nl;
    assert( coords.size() == 2 * num_vertices );
    
    std::vector<std::array<int,3>> tris;
    tris.reserve( num_triangles );
    
    for( int l = Linner; l < Louter; l++ )
    {
        int base_inner  = 3 * ( power_of_two( l-1 ) - 1 ) - 3 * ( power_of_two( Linner-1 ) - 1 );
        int count_inner = 3 * power_of_two( l-1 );
        
        int base_outer  = base_inner + count_inner;
        int count_outer = 2 * count_inner;
        
        for( int i = 0; i < count_inner; i++ )
        {
            
            tris.push_back( { 
                base_inner + i, 
                base_outer + 2*i,
                base_outer + (2*i + 1)%count_outer,
                } );
            tris.push_back( { 
                base_inner + i, 
                base_inner + (i+1)%count_inner,
                base_outer + 2*i + 1,
                } );
            tris.push_back( { 
                base_outer + (2*i+2)%count_outer,
                base_inner + (i+1)%count_inner,
                base_outer + 2*i + 1,
                } );
            
        }
    }
    
    LOG << coords.size() / 2 << space << num_vertices << nl;
    LOG << tris.size() << space << num_triangles << nl;
    assert( tris.size() == num_triangles );
    
//     for( auto t : tris ){ for( auto v : t ) 
//         LOG << v << space; LOG << nl; }
    
    for( auto& t : tris ) //std::sort( t.begin(), t.end() );
    {
      if( t[0] > t[1] ) { int i = t[0]; t[0] = t[1]; t[1] = i; }
      if( t[0] > t[2] ) { int i = t[0]; t[0] = t[2]; t[2] = i; }
      if( t[1] > t[2] ) { int i = t[1]; t[1] = t[2]; t[2] = i; }
      assert( t[0] < t[1] and t[1] < t[2] );
    }
    
//     for( auto t : tris ){ for( auto v : t ) 
//         LOG << v << space; LOG << nl; }
    
    return MeshSimplicial2D(
      2,
      Coordinates( 2, num_vertices, coords ),
      tris
    );
}





inline MeshSimplicial2D SphericalSurface2D( int L = 0 )
{
    
    MeshSimplicial2D ret(
      3,
      Coordinates( 3, 4, {
          std::sqrt(8./9.),               0.0, -1./3., // 0
         -std::sqrt(2./9.),  std::sqrt(2./3.), -1./3., // 1
         -std::sqrt(2./9.), -std::sqrt(2./3.), -1./3., // 2
                       0.0,               0.0,  1.0  , // 3
      } ),
      {
        { 0, 1, 2 },
        { 0, 1, 3 },
        { 0, 2, 3 },
        { 1, 2, 3 }
      }
    );
    
    for( int l = 0; l < L; l++ )
        ret.uniformrefinement();
    
    for( int n = 0; n < ret.getcoordinates().getnumber(); n++ )
    {
        FloatVector point = ret.getcoordinates().getvectorclone( n );
        point.normalize();
        ret.getcoordinates().loadvector( n, point );
    }
    
    return ret;
    
}






#endif
