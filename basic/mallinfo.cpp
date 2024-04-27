#include <stdio.h>

#include "mallinfo.hpp"

#if defined(__linux__) && (__GLIBC__ >= 2) and (__GLIBC_MINOR__ >= 33)

#include <malloc.h>
#include <stdio.h>

void display_mallinfo( bool full )
{
    struct mallinfo2 mi;

    mi = mallinfo2();

    fprintf( stderr, "Heap usage statistics (mallinfo 2):\n");

    if( not full ){
        fprintf( stderr, "Total allocated space (uordblks):      %zu\t\t%.3f GiB\n", (size_t)mi.uordblks, mi.uordblks / 1073741824. );
        fprintf( stderr, "Total free space (fordblks):           %zu\t\t%.3f GiB\n", (size_t)mi.fordblks, mi.fordblks / 1073741824. );
        return;
    }

    fprintf( stderr, "Total non-mmapped bytes (arena):       %zu\t\t%.3f GiB\n", (size_t)mi.arena,    mi.arena    / 1073741824. );
    fprintf( stderr, "# of free chunks (ordblks):            %zu\n",             (size_t)mi.ordblks                             );
    fprintf( stderr, "# of free fastbin blocks (smblks):     %zu\n",             (size_t)mi.smblks                              );
    fprintf( stderr, "# of mapped regions (hblks):           %zu\n",             (size_t)mi.hblks                               );
    fprintf( stderr, "Bytes in mapped regions (hblkhd):      %zu\t\t%.3f GiB\n", (size_t)mi.hblkhd,   mi.hblkhd   / 1073741824. );
    fprintf( stderr, "Max. total allocated space (usmblks):  %zu\t\t%.3f GiB\n", (size_t)mi.usmblks,  mi.usmblks  / 1073741824. );
    fprintf( stderr, "Free bytes held in fastbins (fsmblks): %zu\t\t%.3f GiB\n", (size_t)mi.fsmblks,  mi.fsmblks  / 1073741824. );
    fprintf( stderr, "Total allocated space (uordblks):      %zu\t\t%.3f GiB\n", (size_t)mi.uordblks, mi.uordblks / 1073741824. );
    fprintf( stderr, "Total free space (fordblks):           %zu\t\t%.3f GiB\n", (size_t)mi.fordblks, mi.fordblks / 1073741824. );
    fprintf( stderr, "Topmost releasable block (keepcost):   %zu\t\t%.3f GiB\n", (size_t)mi.keepcost, mi.keepcost / 1073741824. );

    malloc_stats();
    fprintf( stderr, "Glibc version: %d %d\n", __GLIBC__, __GLIBC_MINOR__ );
}

#elif defined(__linux__) && (__GLIBC__ >= 2) and (__GLIBC_MINOR__ >= 15)

#include <malloc.h>
#include <stdio.h>

void display_mallinfo( bool full )
{
    struct mallinfo mi;

    mi = mallinfo();

    fprintf( stderr, "Heap usage statistics (mallinfo 1):\n");

    if( not full ){
        fprintf( stderr, "Total allocated space (uordblks):      %zu\t\t%.3f GiB\n", (size_t)mi.uordblks, mi.uordblks / 1073741824. );
        fprintf( stderr, "Total free space (fordblks):           %zu\t\t%.3f GiB\n", (size_t)mi.fordblks, mi.fordblks / 1073741824. );
        return;
    }

    fprintf( stderr, "Total non-mmapped bytes (arena):       %d\t\t%.3f GiB\n", (int)mi.arena,    mi.arena    / 1073741824. );
    fprintf( stderr, "# of free chunks (ordblks):            %d\n",             (int)mi.ordblks                             );
    fprintf( stderr, "# of free fastbin blocks (smblks):     %d\n",             (int)mi.smblks                              );
    fprintf( stderr, "# of mapped regions (hblks):           %d\n",             (int)mi.hblks                               );
    fprintf( stderr, "Bytes in mapped regions (hblkhd):      %d\t\t%.3f GiB\n", (int)mi.hblkhd,   mi.hblkhd   / 1073741824. );
    fprintf( stderr, "Max. total allocated space (usmblks):  %d\t\t%.3f GiB\n", (int)mi.usmblks,  mi.usmblks  / 1073741824. );
    fprintf( stderr, "Free bytes held in fastbins (fsmblks): %d\t\t%.3f GiB\n", (int)mi.fsmblks,  mi.fsmblks  / 1073741824. );
    fprintf( stderr, "Total allocated space (uordblks):      %d\t\t%.3f GiB\n", (int)mi.uordblks, mi.uordblks / 1073741824. );
    fprintf( stderr, "Total free space (fordblks):           %d\t\t%.3f GiB\n", (int)mi.fordblks, mi.fordblks / 1073741824. );
    fprintf( stderr, "Topmost releasable block (keepcost):   %d\t\t%.3f GiB\n", (int)mi.keepcost, mi.keepcost / 1073741824. );

    malloc_stats();
    fprintf( stderr, "Glibc version: %d %d\n", __GLIBC__, __GLIBC_MINOR__ );
}

#else 

#include <stdio.h>
void display_mallinfo( bool full )
{
    fprintf( stderr, "Heap usage statistics: requested but not available in this implementation.\n");
    #if defined(__GLIBC__) && defined(__GLIBC_MINOR__) 
    fprintf( stderr, "Glibc version: %d %d\n", __GLIBC__, __GLIBC_MINOR__ );
    #else
    fprintf( stderr, "No Glibc version found\n" );
    #endif
}

#endif 
