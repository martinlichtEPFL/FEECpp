#include "logging.hpp"

#ifndef USE_PRIMITIVE_LOGGING
bool log_has_a_fresh_line = true;

// #include <iostream>
// #include <iomanip>
#include <cstdio>
#include <ctime>

// For controlling the floating-point behavior
#include <fenv.h>
#include <float.h>



enum class TextColors : unsigned char
{
    black          = 30,
    black_light    = 90,
    
    red            = 31,
    red_light      = 91,
    
    green          = 32,
    green_light    = 92, 
    
    yellow         = 33, 
    yellow_bright  = 93, 
    
    blue           = 34,
    blue_bright    = 94,
    
    magenta        = 35,
    magenta_bright = 95,
    
    cyan           = 36,
    cyan_light     = 96, 
    
    white          = 37,
    white_light    = 97 
};


Logger::Logger( 
            bool use_cerr, //std::ostream& os,
            const bool do_newline,
            const char* filename,
            const int linenumber
        )
        : 
        use_cerr( use_cerr ), // internalstream( use_cerr ? std::cerr : std::cout ), // internalstream( os ),
        pad_newline_if_there_is_none( do_newline ),
        filename( filename ),
        linenumber( linenumber )
        {
            //*this << std::setprecision(10);
        }


// auto digitalcodenow() -> std::string;

Logger::~Logger()
{
    
    FILE* f = ( use_cerr ? stderr : stdout );
    
    const auto& message_string = internal;
    
    
    // if the log string is empty, just do not print anything 

    if( message_string.empty() ) return;


    // complete the prefix string 

    const std::string time_code = digitalcodenow(); 

    #ifdef USE_COLORED_OUTPUT
    const std::string colorcode_begin = "\033[96m";
    const std::string colorcode_close = "\033[m";
    #else
    const std::string colorcode_begin = "";
    const std::string colorcode_close = "";
    #endif
    const std::string prefix = printf_into_string(
        "%s[%s %s %4d]%s\t", 
        colorcode_begin.c_str(), time_code.c_str(), filename.c_str(), linenumber, colorcode_close.c_str() );
    
    
    // assemble final output string 
    // reserve enough space so that (in principle) each nl can be prefixed

    int number_of_newlines = 0;
    for( char character : message_string ) 
        if( character == '\n' )
            number_of_newlines++;

    std::string output_string;

    output_string.reserve( 1 + message_string.size() + number_of_newlines * prefix.size() );


    // if we have inherited a fresh line, then insert a prefix 
    
    if( log_has_a_fresh_line ){
        log_has_a_fresh_line = false;
        output_string.insert(0,prefix);
    }

    
    // for each nl (except possibly at the final position) insert a prefix after it
    
    for( int c = 0; c < message_string.length(); c++ ) {
    
        char character = message_string[c];

        output_string += character;
        // fputc( character, f );

        if( character == '\n' and c != message_string.length()-1 ) 
            output_string += prefix;
        
    }
    
    
    // if the last character is not a newline and we automatically append a newline,
    // then modify the string accordingly
    if( message_string.back() != nl and pad_newline_if_there_is_none )
        output_string += nl;

    
    log_has_a_fresh_line = ( output_string.back() == nl );
    
    // if( message_string.back() == '\n' or pad_newline_if_there_is_none ) {
    //    
    //     log_has_a_fresh_line = true;
    //  
    //     if( message_string.back() != '\n' and pad_newline_if_there_is_none ) 
    //         output_string += nl;
    //         // fputc( nl, f ); 
    //
    // }

    fputs( output_string.c_str(), f );
    fflush(f);

}






#if defined(_OPENMP)
#include <omp.h>
#endif // #if defined(_OPENMP)

System_Reporter::System_Reporter()
{
    output();
}

System_Reporter::~System_Reporter()
{
    #if defined(_OPENMP)
    // LOG << "###\tSystem Reporter: finished\n";
    #endif
}   

void System_Reporter::output()
{
    
    #if defined(GIT_COMMIT_ID)
    LOG << "###\tCurrent Git commit ID: " << GIT_COMMIT_ID << nl;
    #endif

    time_t t = std::time(nullptr);
    tm tm = *std::localtime(&t);
    LOGPRINTF("###\t%d-%02d-%02d %02d:%02d:%02d\n", tm.tm_year + 1900, tm.tm_mon + 1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
    
    #if defined(_OPENMP)
    LOG << "###\tOpenMP Value: " << _OPENMP << nl;
    LOG << "###\tMaximum number of threads: " << omp_get_max_threads() << nl;
    LOG << "###\tMaximum active levels: " << omp_get_max_active_levels() << nl;
    // LOG << "###\tNested parallelism supported: " << ( omp_get_nested() ? "Yes" : "No" ) << nl;
    LOG << "###\tDynamic adjustment of the number of threads enabled: " << ( omp_get_dynamic() ? "Yes" : "No" ) << nl;
    LOG << "###\tThread limit: " << omp_get_thread_limit() << nl;
    LOG << "###\tNumber of processors available: " << omp_get_num_procs() << nl;
    LOG << "###\tDefault thread affinity: " << ( omp_get_proc_bind() == omp_proc_bind_false ? "No affinity" : "Affinity" ) << nl;
    // LOG << "###\tMaximum number of processors: " << p << " -> " << omp_get_place_num_procs() << nl;
    LOG << "###\tMaximum number of places: " << omp_get_num_places() << nl;
    for( int p = 0; p < omp_get_num_places(); p++ ) {
        LOG << "###\t\tMaximum number of processors per place: " << p << " - > " << omp_get_place_num_procs(p) << nl;
    }
    #else
    LOG << "###\tOpenMP is not enabled.\n";
    #endif

    // #ifdef _OPENMP
    // LOGPRINTF("OpenMP is enabled.\n");
    // LOGPRINTF("Number of threads: %d\n", omp_get_max_threads());
    // LOGPRINTF("Nested parallelism supported: %s\n", omp_get_nested() ? "Yes" : "No");
    // LOGPRINTF("Dynamic adjustment of the number of threads enabled: %s\n", omp_get_dynamic() ? "Yes" : "No");
    // LOGPRINTF("Number of processors available: %d\n", omp_get_num_procs());
    // LOGPRINTF("Default thread affinity: %s\n", omp_get_proc_bind() == omp_proc_bind_false ? "No affinity" : "Affinity");
    // #else
    // LOGPRINTF("OpenMP is not enabled.\n");
    // #endif

    LOGPRINTF("###\tCompiler version: %s\n", __VERSION__);

    #ifdef _DEBUG
    LOGPRINTF("###\tMSVC Debugging flags enabled.\n");
    #else
    LOGPRINTF("###\tMSVC Debugging flags not enabled.\n");
    #endif

    #ifdef NDEBUG
    LOGPRINTF("###\tNDEBUG enabled.\n");
    #else
    LOGPRINTF("###\tNDEBUG not enabled.\n");
    #endif

    #ifdef __OPTIMIZE__
    LOGPRINTF("###\tCompiler optimization level: %d\n", __OPTIMIZE__);
    // Add more performance-related information here
    #endif

    #ifdef __unix__
    LOGPRINTF("###\tOS: Unix\n");
    #elif defined(_WIN32)
    LOGPRINTF("###\tOS: Windows\n");
    #else
    #error "Unknown platform"
    #endif


    #if true and defined(_WIN32)
    LOGPRINTF("###\tFlushing subnormal numbers\n");
    _controlfp_s( nullptr, _DN_FLUSH, _MCW_DN ); // Flush denormals, both operands and results, to zero 
    #endif

}

System_Reporter  omp_reporter;




#endif
