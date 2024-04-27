#ifndef INCLUDEGUARD_EXAMPLES_3D_HPP
#define INCLUDEGUARD_EXAMPLES_3D_HPP


#include "mesh.hpp"
#include "mesh.simplicial3D.hpp"






inline MeshSimplicial3D UnitCube3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 8, {
         0., 0., 0., // 000
         0., 0., 1., // 001
         0., 1., 0., // 010
         0., 1., 1., // 011
         1., 0., 0., // 100
         1., 0., 1., // 101
         1., 1., 0., // 110
         1., 1., 1.  // 111
      } ),
      {
        { 0b000, 0b100, 0b110, 0b111 }, 
        { 0b000, 0b100, 0b101, 0b111 },
        { 0b000, 0b010, 0b110, 0b111 },
        { 0b000, 0b010, 0b011, 0b111 },
        { 0b000, 0b001, 0b101, 0b111 },
        { 0b000, 0b001, 0b011, 0b111 }
      }
    );
}

inline MeshSimplicial3D StandardCube3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 8, {
         -1., -1., -1., // 000
         -1., -1., +1., // 001
         -1., +1., -1., // 010
         -1., +1., +1., // 011
         +1., -1., -1., // 100
         +1., -1., +1., // 101
         +1., +1., -1., // 110
         +1., +1., +1.  // 111
      } ),
      {
        { 0b000, 0b100, 0b110, 0b111 }, 
        { 0b000, 0b100, 0b101, 0b111 },
        { 0b000, 0b010, 0b110, 0b111 },
        { 0b000, 0b010, 0b011, 0b111 },
        { 0b000, 0b001, 0b101, 0b111 },
        { 0b000, 0b001, 0b011, 0b111 }
      }
    );
}



inline MeshSimplicial3D StandardCubeFive3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 8, {
         -1., -1., -1., // 000
         -1., -1., +1., // 001
         -1., +1., -1., // 010
         -1., +1., +1., // 011
         +1., -1., -1., // 100
         +1., -1., +1., // 101
         +1., +1., -1., // 110
         +1., +1., +1.  // 111
      } ),
      {
        { 0b000, 0b110, 0b011, 0b101 }, 
        { 0b000, 0b110, 0b011, 0b010 }, 
        { 0b000, 0b110, 0b100, 0b101 }, 
        { 0b000, 0b001, 0b011, 0b101 },
        { 0b111, 0b110, 0b011, 0b101 } 
      }
    );
}





inline MeshSimplicial3D OctoDiamond3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0., 0., 0., // 0
         1., 0., 0., // 1
         0., 1., 0., // 2
        -1., 0., 0., // 3
         0.,-1., 0., // 4
         0., 0., 1., // 5
         0., 0.,-1.  // 6
      } ),
      {
        { 0, 1, 2, 5 },
        { 0, 2, 3, 5 },
        { 0, 3, 4, 5 },
        { 0, 1, 4, 5 },
        { 0, 1, 2, 6 },
        { 0, 2, 3, 6 },
        { 0, 3, 4, 6 },
        { 0, 1, 4, 6 } 
      }
    );
}




inline MeshSimplicial3D UnitSimplex3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0., 0., 0., // 0
         0., 0., 1., // 1
         0., 1., 0., // 2
         1., 0., 0.  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}

inline MeshSimplicial3D SomeSimplex3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0.2, 0.0,-0.1, // 0
         0.0, 0.3, 0.9, // 1
         0.1, 1.2, 0.1, // 2
         1.5,-0.1,-0.2  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}







inline MeshSimplicial3D RegularSimplex3D()

{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0., 0., 0., // 0
         1., 1., 0., // 1
         1., 0., 1., // 2
         0., 1., 1.  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}

inline MeshSimplicial3D Simplex3D_needle( Float len_x, Float len_y )
{
    assert( len_x >= 0. && len_y >= 0. );
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         0.0,   0.0,   0.0, // 0
         len_x, 0.0,   0.0, // 1
         0.0,   len_y, 0.0, // 2
         0.0,   0.0,   1.0  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}

inline MeshSimplicial3D Simplex3D_tent( Float height )
{
    assert( height >= 0. );
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 4, {
         cos( 0./3. * Constants::twopi ), sin( 0./3. * Constants::twopi ), 0.0, // 0
         cos( 1./3. * Constants::twopi ), sin( 1./3. * Constants::twopi ), 0.0, // 0
         cos( 2./3. * Constants::twopi ), sin( 2./3. * Constants::twopi ), 0.0, // 0
         0.0,   0.0,   height  // 3
      } ),
      {
        { 0, 1, 2, 3 }
      }
    );
}















inline MeshSimplicial3D FicheraCorner3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 17, {
          0.0,  0.0,  0.0, // 0
         
         -1.0, -1.0, -1.0, // 1
          1.0, -1.0, -1.0, // 2
         -1.0,  1.0, -1.0, // 3
          1.0,  1.0, -1.0, // 4
          
         -1.0, -1.0,  1.0, // 5
          1.0, -1.0,  1.0, // 6
         -1.0,  1.0,  1.0, // 7
          // 1.0,  1.0,  1.0, // 8
         
          0.0,  1.0,  1.0, // 8a

         -1.0,  0.0,  0.0, //  9
          1.0,  0.0,  0.0, // 10

          0.0, -1.0,  0.0, // 11
          0.0,  1.0,  0.0, // 12
          
          0.0,  0.0, -1.0, // 13
          0.0,  0.0,  1.0, // 14        

          1.0,  0.0,  1.0, // 8b // 15
          1.0,  1.0,  0.0, // 8c // 16
      } ),
      {  
        { 0, 1, 3,  9 }, { 0, 3, 7,  9 }, { 0, 7, 5,  9 }, { 0, 5, 1,  9 },
        { 0, 2, 4, 10 }, { 0, 4,16, 10 }, { 0,15, 6, 10 }, { 0, 6, 2, 10 }, // 8c, 8b
        { 0, 1, 2, 11 }, { 0, 2, 6, 11 }, { 0, 6, 5, 11 }, { 0, 5, 1, 11 },
        { 0, 3, 4, 12 }, { 0, 4,16, 12 }, { 0, 8, 7, 12 }, { 0, 7, 3, 12 }, // 8c, 8(a)
        { 0, 1, 2, 13 }, { 0, 2, 4, 13 }, { 0, 4, 3, 13 }, { 0, 3, 1, 13 },
        { 0, 5, 6, 14 }, { 0, 6,15, 14 }, { 0, 8, 7, 14 }, { 0, 7, 5, 14 }  // 8b, 8(a)
      }
    );
}





inline MeshSimplicial3D CrissCrossCube3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 15, {
          0.0,  0.0,  0.0, // 0
         
         -1.0, -1.0, -1.0, // 1
          1.0, -1.0, -1.0, // 2
         -1.0,  1.0, -1.0, // 3
          1.0,  1.0, -1.0, // 4
          
         -1.0, -1.0,  1.0, // 5
          1.0, -1.0,  1.0, // 6
         -1.0,  1.0,  1.0, // 7
          1.0,  1.0,  1.0, // 8
         
         -1.0,  0.0,  0.0, //  9
          1.0,  0.0,  0.0, // 10

          0.0, -1.0,  0.0, // 11
          0.0,  1.0,  0.0, // 12
          
          0.0,  0.0, -1.0, // 13
          0.0,  0.0,  1.0, // 14        
      } ),
      {
        { 0, 1, 3,  9 }, { 0, 3, 7,  9 }, { 0, 7, 5,  9 }, { 0, 5, 1,  9 },
        { 0, 2, 4, 10 }, { 0, 4, 8, 10 }, { 0, 8, 6, 10 }, { 0, 6, 2, 10 },
        { 0, 1, 2, 11 }, { 0, 2, 6, 11 }, { 0, 6, 5, 11 }, { 0, 5, 1, 11 },
        { 0, 3, 4, 12 }, { 0, 4, 8, 12 }, { 0, 8, 7, 12 }, { 0, 7, 3, 12 },
        { 0, 1, 2, 13 }, { 0, 2, 4, 13 }, { 0, 4, 3, 13 }, { 0, 3, 1, 13 },
        { 0, 5, 6, 14 }, { 0, 6, 8, 14 }, { 0, 8, 7, 14 }, { 0, 7, 5, 14 } 
      }
    );
}




inline MeshSimplicial3D CrossedBricks3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 20, {
         0.0,  0.0,  0.0, // 0
         
         // x    y      z
         
         0., 0., 1., // 1           *
         0., 1., 0., // 2
         0., 1., 1., // 3
         
         1., 0., 0., // 4
         1., 0., 1., // 5
         1., 1., 0., // 6
         1., 1., 1.,// 7

        -1., 0., 0., // 8
        -1., 0., 1., // 9           * 
        -1., 1., 0., // 10
        -1., 1., 1., // 11

         0., -1., 0., // 12 <- 2 
         0., -1., 1., // 13 <- 3    *
        -1., -1., 0., // 14 <- 10
        -1., -1., 1., // 15 <- 11   *

         0.,  0., -1., // 16 <- 1  
        -1.,  0., -1., // 17 <- 9  
         0., -1., -1., // 18 <- 13 
        -1., -1., -1.  // 19 <- 15 



      } ),
      {
        { 0,  4,  6,  7 }, 
        { 0,  4,  5,  7 },
        { 0,  2,  6,  7 },
        { 0,  2,  3,  7 },
        { 0,  1,  5,  7 },
        { 0,  1,  3,  7 },
      
        { 0,  8, 10, 11 }, 
        { 0,  8,  9, 11 },
        { 0,  2, 10, 11 },
        { 0,  2,  3, 11 },
        { 0,  1,  9, 11 },
        { 0,  1,  3, 11 },
      
        { 0,  8, 14, 15 }, 
        { 0,  8,  9, 15 },
        { 0, 12, 14, 15 },
        { 0, 12, 13, 15 },
        { 0,  1,  9, 15 },
        { 0,  1, 13, 15 },
      
        { 0,  8, 14, 19 }, 
        { 0,  8, 17, 19 },
        { 0, 12, 14, 19 },
        { 0, 12, 18, 19 },
        { 0, 16, 17, 19 },
        { 0, 16, 18, 19 }
      }
    );
}




inline MeshSimplicial3D CrossedBricks_Five3D()
{
    return MeshSimplicial3D(
      3,  
      Coordinates( 3, 20, {
         -1.,  0.,  1., // 0000
         -1.,  0.,  0., // 0001
         -1., +1.,  1., // 0010
         -1., +1.,  0., // 0011
          0.,  0.,  1., // 0100
          0.,  0.,  0., // 0101
          0., +1.,  1., // 0110
          0., +1.,  0., // 0111

         +1.,  0.,  1., // 1000
         +1.,  0.,  0., // 1001
         +1., +1.,  1., // 1010
         +1., +1.,  0., // 1011

         -1., -1.,  1., // 1100
         -1., -1.,  0., // 1101
          0., -1.,  1., // 1110
          0., -1.,  0., // 1111

         -1.,  0., -1., // 16 <- 0000
          0.,  0., -1., // 17 <- 0100
         -1., -1., -1., // 18 <- 1100
          0., -1., -1.  // 19 <- 1110
      } ),
      {
        { 0b0000, 0b0110, 0b0011, 0b0101 }, 
        { 0b0000, 0b0110, 0b0011, 0b0010 }, 
        { 0b0000, 0b0110, 0b0100, 0b0101 }, 
        { 0b0000, 0b0001, 0b0011, 0b0101 },
        { 0b0111, 0b0110, 0b0011, 0b0101 }, 

        { 0b1000, 0b0110, 0b1011, 0b0101 }, 
        { 0b1000, 0b0110, 0b1011, 0b1010 }, 
        { 0b1000, 0b0110, 0b0100, 0b0101 }, 
        { 0b1000, 0b1001, 0b1011, 0b0101 },
        { 0b0111, 0b0110, 0b1011, 0b0101 }, 

        { 0b0000, 0b1110, 0b1101, 0b0101 }, 
        { 0b0000, 0b1110, 0b1101, 0b1100 }, 
        { 0b0000, 0b1110, 0b0100, 0b0101 }, 
        { 0b0000, 0b0001, 0b1101, 0b0101 },
        { 0b1111, 0b1110, 0b1101, 0b0101 }, 

        {     16,     19, 0b1101, 0b0101 }, 
        {     16,     19, 0b1101,     18 }, 
        {     16,     19,     17, 0b0101 }, 
        {     16, 0b0001, 0b1101, 0b0101 },
        { 0b1111,     19, 0b1101, 0b0101 }  
      }
    );
}







#endif