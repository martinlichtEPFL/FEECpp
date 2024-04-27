#ifndef INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP
#define INCLUDEGUARD_UTILITY_CONVERGENCETABLE_HPP

// #include <ostream>
#include <string>
#include <vector>

#include "../basic.hpp"


class ConvergenceTable
{
    typedef long double EntryType; 
    
    private:
        
        std::vector<std::vector<EntryType>> entries;
        std::vector<std::string> seriesheaders;
        bool make_new_row;
        
    public:
        
        std::string table_name;
        bool display_convergence_rates;
        bool print_rowwise_instead_of_columnwise;
    
    public:

        explicit ConvergenceTable( const std::string& table_name = "---------- Default Table Name ----------" );

        ConvergenceTable( const ConvergenceTable& ) = default;
        ConvergenceTable( ConvergenceTable&& ) = default;
        ConvergenceTable& operator=( const ConvergenceTable& ) = default;
        ConvergenceTable& operator=( ConvergenceTable&& ) = default;
        
        
        void insert_numerical_entry( EntryType entry );
        void insert_seriesheader( const std::string& seriesheader );
        void insert_newline();

        Float get_convergence_rate( int row_index, int column_index );


        std::string text() const;
        
        void lg() const;
        
        void lg( bool display_convergence_rates ) const;
        
//         void print( std::ostream& os );
        

        std::string text( bool display_convergence_rates ) const;
        
        std::string text_standard( bool display_convergence_rates ) const;

        std::string text_transpose( bool display_convergence_rates ) const;


        std::string TeXtabular() const;
        std::string TeXtabular( const std::vector<bool>& ) const;


        // stream operators 

        ConvergenceTable& operator<<( float entry );

        ConvergenceTable& operator<<( double entry );

        ConvergenceTable& operator<<( long double entry );

        ConvergenceTable& operator<<( const std::string& seriesheader );

        ConvergenceTable& operator<<( char code );

};




        


#endif
