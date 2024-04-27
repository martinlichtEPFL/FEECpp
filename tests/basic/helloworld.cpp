
#include <iostream>
#include <new>

#include "../../basic/mallinfo.hpp"

using namespace std;

/*
The purpose of this program is to check whether your development environment is set up correctly,
including compiler, linker, debug software, etc...
It also outputs the C++ version.
*/

int main( int argc, char *argv[] )
{
    cout << "Hello World! " << endl;
                
    // C++ standard
    cout << "C++ version: " << __cplusplus << endl;
        
    // Compiler identification
    #if defined(__clang__)
    cout << "Compiler: Clang " << __clang_version__ << endl;
    #elif defined(__GNUC__)
    cout << "Compiler: GCC " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << endl;
    #elif defined(_MSC_VER)
    cout << "Compiler: MSVC " << _MSC_VER << endl;
    #else
    cout << "Compiler: unknown" << endl;
    #endif
    
    // Standard library and relevant macros
    #if defined(_LIBCPP_VERSION)
    cout << "Standard Library: LLVM libc++ " << _LIBCPP_VERSION << endl;
    #elif defined(__GLIBCXX__)
    cout << "Standard Library: GNU libstdc++ " << __GLIBCXX__ << endl;
    #elif defined(_MSVC_STL_VERSION)
    cout << "Standard Library: MSVC STL " << _MSVC_STL_VERSION << endl;
    #endif
    
    // GNU libc version (for GNU/Linux systems)
    #ifdef __GLIBC__
    cout << "GNU libc version: " << __GLIBC__ << "." << __GLIBC_MINOR__ << endl;
    #endif

    // Operating system identification
    #if defined(__linux__)
    cout << "Operating system: Linux" << endl;
    #elif defined(__APPLE__)
    cout << "Operating system: Apple" << endl;
    #elif defined(__FreeBSD__) or defined(__NetBSD__) or defined(__OpenBSD__)
    cout << "Operating system: BSD" << endl;
    #elif defined(_WIN32) or defined(_WIN64)
    cout << "Operating system: Windows" << endl;
    #else
    cout << "Operating system: unknown" << endl;
    #endif
    
    display_mallinfo();
    
    cout << "Now an intentional leak..." << endl;
    
    int * p = new (nothrow) int[10000];
    
    cout << p[8] << endl;
    
    return 0;
}

