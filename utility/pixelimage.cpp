
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <string>


#include "../external/stb_image.h"
#include "../external/stb_image_write.h"

#include "pixelimage.hpp"



// // utility for getting file ending
// std::string tail(std::string const& source, size_t const length) {
//   if( length >= source.size()) { return source; }
//   return source.substr(source.size() - length);
// } 




PixelImage readPixelImage( const std::string& str )
{

    int width, height, channels;
    
    LOG << "Loading: " << str << nl;

    unsigned char *image = stbi_load( str.c_str(), &width, &height, &channels, 0 );
    
    if( image == nullptr ) { 
        LOG << "STATUS: " << stbi_failure_reason() << nl;
        Assert( image != nullptr );
    }
    
    // std::cout << "Loaded image with a width of " << width << ", a height of " << height << " and " << channels << " channels" << std::endl;

    PixelImage ret( height, width ); 
    
    for( int row = 0; row < height; row++) 
    for( int col = 0; col < width; col++)
    {
        ret(row,col).red   = image[ channels * ( row * width + col ) + 0 ];
        ret(row,col).green = image[ channels * ( row * width + col ) + 1 ];
        ret(row,col).blue  = image[ channels * ( row * width + col ) + 2 ];
    }

    stbi_image_free( image );

    return ret;

}

void savePixelImage( const PixelImage& pim, const std::string& str )
{
    
    const int width    = pim.getwidth();
    const int height   = pim.getheight();
    const int channels = 3;

    unsigned char *image = new unsigned char[ width * height * channels ];

    for( int row = 0; row < height; row++ ) 
    for( int col = 0; col < width;  col++ )
    {
        image[ channels * ( row * width + col ) + 0 ] = pim(row,col).red;
        image[ channels * ( row * width + col ) + 1 ] = pim(row,col).green;
        image[ channels * ( row * width + col ) + 2 ] = pim(row,col).blue;
    }

    const bool output_success = stbi_write_png( str.c_str(), width, height, channels, image, width * channels );
    
    Assert( output_success );
    
    delete[] image;
    
}








    // // Initialize ImageMagick library
    // InitializeMagick(nullptr);

    // // Load the image file
    // Image image;
    
    // image.read("lena_color.tiff");

    // // Get the dimensions of the image
    // size_t width = image.columns();
    // size_t height = image.rows();

    // // Allocate stack-allocated arrays to store the pixel values
    // unsigned char pixels_r[height][width];
    // unsigned char pixels_g[height][width];
    // unsigned char pixels_b[height][width];

    // // Get the pixel values and store them in the arrays
    // Pixels view(image);
    // PixelPacket* pixel_cache = view.get(0, 0, width, height);
    // for( size_t row = 0; row < height; row++) {
    //     for( size_t col = 0; col < width; col++) {
    //         size_t index = row * width + col;
    //         pixels_r[row][col] = pixel_cache[index].red;
    //         pixels_g[row][col] = pixel_cache[index].green;
    //         pixels_b[row][col] = pixel_cache[index].blue;

    //         pixel_cache[index].red   = 255 - pixel_cache[index].red;
    //         pixel_cache[index].green = 255 - pixel_cache[index].green;
    //         pixel_cache[index].blue  = 255 - pixel_cache[index].blue;
    //     }
    // }

    // view.sync();

    // image.write("enal_color.tiff");

    // // Print the pixel values
    // for( size_t row = 0; row < height; row++) {
    //     for( size_t col = 0; col < width; col++) {
    //         std::cout << "(" << (int)pixels_r[row][col] << "," << (int)pixels_g[row][col] << "," << (int)pixels_b[row][col] << ") ";
    //     }
    //     std::cout << std::endl;
    // }

    // // Terminate ImageMagick library    
    // return 0;
// }


// RGBstring rgb_to_string( unsigned char r, unsigned char g, unsigned char b )
// {
//     char result[8] = {0,0,0,0,0,0,0,0};
//     sprintf(result, "#%02x%02x%02x", r, g, b);
//     return std::string( result );
// }

// unsigned char[3] string_to_rgb( RGBstring rgb )
// {
//     unsigned char r, g, b;
//     r = ( rgb[1] - 'A' ) * 16 + ( rgb[2] - 'A' );
//     g = ( rgb[3] - 'A' ) * 16 + ( rgb[4] - 'A' );
//     b = ( rgb[5] - 'A' ) * 16 + ( rgb[6] - 'A' );
//     return { r, g, b };
// }

// int leading_digits( double num )
// {
//     // If the number is negative, take its absolute value
//     // Calculate the order of magnitude of the number
//     // If the number is greater than 1, return the integer part of the order, else, return 0
    
//     if( num < 0 ) num = -num;
//     int order = (int)ceil( log10(num) );
//     if( order >= 0 ) return order + 1; else return 0;
// }

// std::string render_number( double num, int tail )
// {
//     int lead = leading_digits( num );
//     char str[1+lead+1+tail+1];
//     sprintf( str, "%0*.*f", lead+1+tail, tail, num);
//     return std::string(str);
// }


/*
    <polygon points="10,20 30,40 50,60" fill="#FF3399" stroke="none"></polygon>
    <polygon points="10,20 40,50 0010,10" fill="#FF3399" stroke="#11CCDD" stroke-width="1"></polygon>

    #include <iostream>
    #include <sstream>
    #include <string>

    std::string generatePolygonString(
        double x1, double y1, double x2, double y2, double x3, double y3, 
        std::string fill, 
        std::string stroke,
        double stroke_width
    ) {
        std::stringstream ss;
        ss << "<polygon points=\"" 
           << x1 << "," << y1 << " " 
           << x2 << "," << y2 << " " 
           << x3 << "," << y3 << "\""  
        ss << " fill=\"" << fill << "\""  
        ss << " stroke=\"" << stroke << "\""  
        ss << " stroke-width=\"" << stroke_width << "\""  
        ss << "></polygon>";
        return ss.str();
    }

    int main() {
        std::string polygonString = generatePolygonString(10, 20, 30, 40, 50, 60, "#FF3399", "none");
        std::cout << polygonString << std::endl;
        return 0;
    }


    void writeSVG( 
        const std::string filename,
        const std::vector<std::array<3,int>>& triangles, 
        const Coordinates& coords,
        const FloatVector* triangle_values = nullptr
    ) {



    };

    int main() {
        unsigned char r = 255;
        unsigned char g = 0;
        unsigned char b = 128;
        RGBstring rgb = rgb_string( r, g, b );
        printf( "RGB string: %s\n", rgb );
        return 0;
    }

*/

