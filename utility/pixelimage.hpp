#ifndef INCLUDEGUARD_UTILITY_PIXELIMAGE_HPP
#define INCLUDEGUARD_UTILITY_PIXELIMAGE_HPP

#include <string>

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <vector>

#include "../basic.hpp"

enum class ColorChannel : unsigned char { red = 0, green = 1, blue = 2 };

struct PixelColor
{
    unsigned char red;
    unsigned char green;
    unsigned char blue;

    unsigned char operator()(ColorChannel cc) const { 
        switch(cc) {
            case ColorChannel::red:   return red;
            case ColorChannel::green: return green;
            case ColorChannel::blue:  return blue;
            default: unreachable();
        }
        unreachable();
    }
};

class PixelImage
{
    private:
        size_t height;
        size_t width;
        std::vector<PixelColor> data;

    public:

        PixelImage( size_t height, size_t width ) 
        : height(height), width(width), data( width*height, {128,128,128} ) 
        {}
        
        PixelImage( const PixelImage& ) = default;
        PixelImage( PixelImage&& ) = default;
        
        ~PixelImage() {}

        size_t getheight() const { return height; }

        size_t getwidth() const { return width; }

        PixelColor& operator()( size_t row, size_t col ) {
            assert( row >= 0 && row < height && col >= 0 && col < width );
            return data.at( row * width + col );
        }
        
        PixelColor operator()( size_t row, size_t col ) const {
            assert( row >= 0 && row < height && col >= 0 && col < width );
            return data.at( row * width + col );
        }

        auto get_interpolated_function() const {
            return [&]( double x, double y, ColorChannel cc ) -> double
            {
                const PixelImage& pixelimage = *this;
                
                assert( 0 <= y and y < height );
                assert( 0 <= x and x < width  );
                
                int ny = static_cast<int>(floor( y ));
                int nx = static_cast<int>(floor( x ));

                const double lx = x - nx;
                const double ly = y - ny;

                assert( 0. <= ly and ly <= 1. );
                assert( 0. <= lx and lx <= 1. );
                assert( 0 <= ny and ny < height );
                assert( 0 <= nx and nx < width  );
                // assert( 0 <= ny and ny+1 < height );
                // assert( 0 <= nx and nx+1 < width  );

                // LOG << nx << space << ny << space << width << space << height << nl;

                PixelColor p00 = pixelimage( ny+0, nx+0 );
                // PixelColor p01 = pixelimage( ny+0, nx+1 );
                // PixelColor p10 = pixelimage( ny+1, nx+0 );
                // PixelColor p11 = pixelimage( ny+1, nx+1 );

                return p00(cc);
                // return       ly *     lx * p00(cc) 
                //        + (1-ly) *     lx * p01(cc) 
                //        +     ly * (1-lx) * p10(cc) 
                //        + (1-ly) * (1-lx) * p11(cc);
            };
        }

        auto get_interpolated_red() const {
            auto rgb = get_interpolated_function();
            return [rgb]( double x, double y ) -> double { return rgb(x,y,ColorChannel::red); };
        }
        
        auto get_interpolated_green() const {
            auto rgb = get_interpolated_function();
            return [rgb]( double x, double y ) -> double { return rgb(x,y,ColorChannel::green); };
        }

        auto get_interpolated_blue() const {
            auto rgb = get_interpolated_function();
            return [rgb]( double x, double y ) -> double { return rgb(x,y,ColorChannel::blue); };
        }
        
        
};


PixelImage readPixelImage( const std::string& str ); 

void savePixelImage( const PixelImage& pim, const std::string& str );




#endif
