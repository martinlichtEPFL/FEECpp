#ifndef INCLUDEGUARD_LOGGING_HPP
#define INCLUDEGUARD_LOGGING_HPP

#include "basic.hpp"



#ifndef USE_PRIMITIVE_LOGGING

#include <cstdio>

// #include <ostream>
#include <string>
// #include <sstream>

// #include "logger.hpp"
// #include "prefixbuffer.hpp"

extern bool log_has_a_fresh_line;


// This variable has an instance in every translation unit 
// It is not global for the entire program 
class Logger //: public std::ostringstream
{
    private:
        
        std::string internal = "";
        
        bool use_cerr; //std::ostream& internalstream;
        bool pad_newline_if_there_is_none;
        std::string filename;
        int linenumber;
    
        // static bool log_has_a_fresh_line = true;

        // static const bool print_file_and_line = false;
        
    public:
    
        explicit 
        // inline
        Logger( 
            bool use_cerr, //std::ostream& os,
            const bool do_newline = false,
            const char* filename = "UNKNOWN",
            const int linenumber = -1
        )
        ;
        // : 
        // internalstream( os ),
        // pad_newline_if_there_is_none( do_newline ),
        // filename( filename ),
        // linenumber( linenumber )
        // {}


        Logger& operator<<( char input ) {
            internal += input;
            return *this;
        }

        Logger& operator<<( const char* input ) {
            internal += input;
            return *this;
        }

        Logger& operator<<( const std::string& input ) {
            internal += input;
            return *this;
        }

        Logger& operator<<( const std::string&& input ) {
            internal += input;
            return *this;
        }

        Logger& operator<<( const void* input ) {
            char buffer[ sizeof(decltype(input)) * 2 + 10 + 1 ]; // how pointers are printed is implementation-defined 
            std::snprintf( buffer, sizeof(buffer), "%p", input );
            internal += buffer;
            return *this;
        }

        template <typename T>
        typename std::enable_if< std::is_integral<T>::value && std::is_signed<T>::value, Logger&>::type
        operator<<(T input) {
            char buffer[ std::numeric_limits<T>::digits10+1 + 1 + 1];
            std::snprintf(buffer, sizeof(buffer), "%jd", static_cast<intmax_t>(input));
            internal += buffer;
            return *this;
        }
    
        template <typename T>
        typename std::enable_if< std::is_integral<T>::value && std::is_unsigned<T>::value, Logger&>::type
        operator<<(T input) {
            char buffer[ std::numeric_limits<T>::digits10+1 + 1 + 1];
            std::snprintf(buffer, sizeof(buffer), "%ju", static_cast<uintmax_t>(input));
            internal += buffer;
            return *this;
        }
    
        // template <typename T, typename = decltype(std::declval<T>().text())>
        // Logger& operator<<(const T& input) {
        //     std::string text = input.text();
        //     internal += text;
        //     return *this;
        // }
    
        // Logger& operator<<( intmax_t input ) {
        //     char buffer[ sizeof(decltype(input)) * 3 + 1 + 1 ];
        //     std::snprintf( buffer, sizeof(buffer), "%jd", input );
        //     internal += buffer;
        //     return *this;
        // }

        // Logger& operator<<( uintmax_t input ) {
        //     char buffer[ sizeof(decltype(input)) * 3 + 1 + 1 ];
        //     std::snprintf( buffer, sizeof(buffer), "%ju", input );
        //     internal += buffer;
        //     return *this;
        // }

        template <typename T>
        typename std::enable_if<std::is_floating_point<T>::value, Logger&>::type
        operator<<( T input ) {
            char buffer[ std::numeric_limits<float>::max_digits10 + std::numeric_limits<float>::max_exponent10 + 10 + 1];
            std::snprintf( buffer, sizeof(buffer), "%.10Lg", (long double)input );
            internal += buffer;
            return *this;
        }

        ~Logger();

};






// class Logger2
// {
//     private:
//         std::ostream& internalstream;
//         std::string prefix;
//         bool pad_newline_if_there_is_none;
//         std::string filename;
//         int linenumber;
    
//         bool print_file_and_line = false;
//         std::ostringstream internalbuffer;
        
//     public:
    
//         explicit inline Logger2( 
//             std::ostream& os,
//             const std::string& prefix = "",
//             const bool do_newline = false,
//             const char* filename = "UNKNOWN",
//             const int linenumber = -1
//         )
//         : internalstream( os ), prefix( prefix ), pad_newline_if_there_is_none( do_newline )
//         {}

//         inline ~Logger2()
//         {
            
//             const auto str = internalbuffer.str();
//             ................. internalstream << str;
            
//         }

//         template<class T>
//         Logger& operator<<( const T& t )
//         {
//             internalbuffer << t;
//             return *this;
//         }
        
//         Logger& operator<<( std::ostream& (*const f)(std::ostream&) )
//         {
//             f( internalbuffer );
//             return *this;
//         }
        
// };





// template< typename L, typename... Params >
// void printf_into_logger( L logger, const char* formatstring, Params... args )
// {
//     logger << printf_into_string( formatstring, args... );
    
//     // std::size_t length = std::snprintf(nullptr, 0, formatstring, args... ) + 1;
//     // char* str = new char[length];
//     // std::snprintf( str, length, formatstring, args... );
    
//     // logger << str;

//     // delete[] str;
// }

#endif // USE_PRIMITIVE_LOGGING






#ifndef USE_PRIMITIVE_LOGGING
// returns a temporary logger to write stuff to, and line breaks on destruction 
// Example usage:
//     LOG << "This is a short message with a number: " << 5;      
//     ERR << "This is an error message.";      

#define LOG     Logger( false, false, __FILE__, __LINE__ )
#define ERR     Logger( true, false, __FILE__, __LINE__ )

#else 

#include <iostream>

#define LOG     std::cout
#define ERR     std::cerr

#endif // USE_PRIMITIVE_LOGGING



// utilize the printf template for stream-like objects 


#if __cplusplus < 202002L
#define LOGPRINTF(...) printf_into_stream( LOG, __VA_ARGS__ );
#define ERRPRINTF(...) printf_into_stream( ERR, __VA_ARGS__ );
#else
#define LOGPRINTF( formatstring, ...) printf_into_stream( LOG, formatstring __VA_OPT__(,) __VA_ARGS__ );
#define ERRPRINTF( formatstring, ...) printf_into_stream( ERR, formatstring __VA_OPT__(,) __VA_ARGS__ );
#endif



// treat the following macros as PRINT 'str' commands
// Example usage:
//     NOTE "This is a note"
//     WARN "This is a warning"
//     ALERT "This is an alert"
//     ERROR "This is an error"

#ifndef USE_PRIMITIVE_LOGGING

#define NOTE    Logger( false, true, __FILE__, __LINE__ ) <<
#define NOTICE  Logger( false, true, __FILE__, __LINE__ ) <<

#define WARNING Logger( true, true, __FILE__, __LINE__ ) <<
#define ALERT   Logger( true, true, __FILE__, __LINE__ ) <<
#define ERROR   Logger( true, true, __FILE__, __LINE__ ) <<

#else 

#define NOTE    std::cout <<
#define NOTICE  std::cout <<

#define WARNING std::cerr <<
#define ALERT   std::cerr <<
#define ERROR   std::cerr <<

#endif // USE_PRIMITIVE_LOGGING


// emit the current file and line number into the log stream 
// Example usage:
//     PING;

#define PING LOG << "PING: " << __FILE__ << ":" << __LINE__ << nl;







////////////////////////////////////////////
// 
//      logging via variadic templates
// 
////////////////////////////////////////////

// inline void lg(){}
// 
// // template<typename T>
// // inline void lg( T arg )
// // {
// //     LOG << arg << nl;
// // }
// 
// template<typename T, typename... Ts>
// inline void lg( const T arg, const Ts... args )
// {
//     LOG << arg << nl;
//     lg( args... );
// }







////////////////////////////////////////////
// 
//      OMP reporter
// 
////////////////////////////////////////////

// #ifdef _OPENMP

struct System_Reporter
{
    System_Reporter();
    ~System_Reporter();
    void output();
};

extern System_Reporter  omp_reporter;

// #endif // _OPENMP








// LEGACY DEFINITIONS:

// inline void ping() { std::clog << "ping" << std::endl; }
// inline void pong() { std::clog << "pong" << std::endl; }
// inline void peng() { std::clog << "peng" << std::endl; }
// inline void pang() { std::clog << "pang" << std::endl; }
// inline void pung() { std::clog << "pung" << std::endl; }
// 
// 
// static std::ostream* lognotice = &std::clog;
// static std::ostream* loginfo   = &std::clog;
// 
// static std::ostream* logwarn   = &std::cerr;
// static std::ostream* logalert  = &std::cerr;
// static std::ostream* logerr    = &std::cerr;

#endif // INCLUDEGUARD_LOGGING_HPP
