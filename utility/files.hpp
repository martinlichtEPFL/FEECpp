#ifndef INCLUDEGUARD_UTILITY_FILES_HPP
#define INCLUDEGUARD_UTILITY_FILES_HPP


#include <fstream>
#include <string>

#include "../basic.hpp"





inline std::string get_parent_directory( const std::string& filepath ) 
{
    size_t pos = filepath.find_last_of("/\\");
    
    if (pos == std::string::npos) {
        // No directory separator found, implying the file is in the current directory
        return ".";
    } else if (pos == 0) {
        // The file is in the root directory
        return "/";
    } else {
        // Extract the parent directory path
        return filepath.substr(0, pos);
    }
}



inline std::fstream openinputfile( const std::string& filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::in );
    Assert( myfile.is_open(), filename );
    return myfile;
}

inline std::fstream openoutputfile( const std::string& filename )
{
    std::fstream myfile;
    myfile.open(filename, std::ios::out );
    Assert( myfile.is_open(), filename );
    return myfile;
}

inline bool fileexists( const std::string& filename )
{
    std::ifstream file( filename.c_str() );
    return file.good();
}

inline std::string adaptfilename( std::string filename )
{
    if( !fileexists( filename ) ) 
        return filename;
    
    auto dotpos = filename.rfind( "." );
    std::string pre_dot  = filename.substr(0,dotpos);
    std::string post_dot = filename.substr(dotpos+1);
    
    for( int i = 0; fileexists( filename = pre_dot + std::string(".") + std::to_string(i) + std::string(".") + post_dot ); i++ );
    
    return filename;
}


inline std::string experimentfile( const std::string& basename, const std::string& extension = "vtk" )
{
    std::string ret;
    for( int i = 0; fileexists( ret = basename + std::string(".") + std::to_string(i) + std::string(".") + std::string(extension) ); i++ );
    return ret;
}

inline std::string getbasename( const std::string& path )
{
      std::string::size_type last_slash = path.find_last_of("/");
      std::string::size_type last_dot   = path.find_last_of(".");
      
      Assert( last_dot != std::string::npos );
      
      std::string::size_type begin = ( last_slash != std::string::npos ) ? ( last_slash+1 ) : 0;
      std::string::size_type end   = last_dot;
      
      return path.substr(begin,end-begin);
}



#endif
