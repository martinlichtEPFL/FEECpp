
#include <cstdlib>
#include <cstdio>
#include <cstddef>
#include <cstdint>

void display_mallinfo();

#if 0 
void* operator new( decltype(sizeof(0)) size )
{
    if( size == 0 )
    {
        ++size; // avoid std::malloc(0) which may return nullpointer on success
        std::puts("Called for allocation of zero bytes");
    }
 
    if( size < 0 )
    {
        std::fprintf( stderr, "\n\nAllocating negative memory size: %zu\n\n", size );
    }

    if( size >= 536870912 ) // half a gigabyte 
    // if( size >= 1024  ) // kilobyte  
    {
        std::fprintf( stderr, "\n\nAllocating large chunk: %zu\n\n", size );
        std::fprintf( stderr, "SIZE_MAX = %zu\t\t%.3f GiB\n\n", SIZE_MAX, SIZE_MAX / 1073741824. );
    }
 
    void *pointer = std::malloc(size);
    if( pointer )
        return pointer;
 
    display_mallinfo();

    std::fprintf( stderr, "MEMORY ALLOCATION FAILED\n" );
    std::fprintf( stderr, "\t size = %zu\n", size );
    std::fprintf( stderr, "\t pointer = %p", pointer );
    std::fprintf( stderr, "GOOD BYE!\n" );
    exit(1);
}
 
// void* operator new( decltype(sizeof(0)) size, decltype(sizeof(0)) align )
// {
//     if( size == 0 )
//     {
//         ++size; // avoid std::malloc(0) which may return nullpointer on success
//         std::puts("Called for allocation of zero bytes");
//     }
 
//     if( size < 0 )
//     {
//         std::fprintf( stderr, "\n\nAllocating negative memory size: %zu\n\n", size );
//     }

//     if( size >= 536870912 ) // half a gigabyte 
//     // if( size >= 1024  ) // kilobyte  
//     {
//         std::fprintf( stderr, "\n\nAllocating large chunk: %zu\n\n", size );
//         std::fprintf( stderr, "SIZE_MAX = %zu\t\t%.3f GiB\n\n", SIZE_MAX, SIZE_MAX / 1073741824. );
//         std::fprintf( stderr, "\t align = %zu\n", align );
//     }
 
//     void *pointer = std::malloc(size,align);
//     if( pointer )
//         return pointer;
 
//     display_mallinfo();

//     std::fprintf( stderr, "MEMORY ALLOCATION FAILED\n" );
//     std::fprintf( stderr, "\t size = %zu\n", size );
//     std::fprintf( stderr, "\t align = %zu\n", align );
//     std::fprintf( stderr, "\t pointer = %p", pointer );
//     std::fprintf( stderr, "GOOD BYE!\n" );
//     exit(1);
// }
 
void operator delete(void* pointer) noexcept
{
    // std::puts("global op delete called\n");
    if( pointer == nullptr ) return;
    std::free(pointer);
}

#endif

