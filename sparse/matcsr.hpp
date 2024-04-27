#ifndef INCLUDEGUARD_SPARSE_MATCSR_HPP
#define INCLUDEGUARD_SPARSE_MATCSR_HPP

#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "sparsematrix.hpp"




/************************
****
****  Class for Sparse Matrices in CSR format
****  - instantiates LinearOperator
****  
************************/




class MatrixCSR:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        /* Constructors */
        
        explicit MatrixCSR( int rows, int columns, 
                            const std::vector<int>& A, 
                            const std::vector<int>& C, 
                            const std::vector<Float>& V );

        explicit MatrixCSR( int rows, int columns, 
                            const std::vector<int>&& A, 
                            const std::vector<int>&& C, 
                            const std::vector<Float>&& V );

        explicit MatrixCSR( const SparseMatrix& mat );

        explicit MatrixCSR( int rows, int columns );

        /* standard interface */ 
        
        MatrixCSR() = delete;
        MatrixCSR( const MatrixCSR& );
        MatrixCSR& operator=( const MatrixCSR& );
        MatrixCSR( MatrixCSR&& );
        MatrixCSR& operator=( MatrixCSR&& );
        virtual ~MatrixCSR( );

        
        /* standard methods for operators */
        
        virtual void check() const override;
        virtual std::string text() const override;
        // virtual void printplain( std::ostream& ) const;


        /* OTHER METHODS */
        
        virtual MatrixCSR* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override;
        
        
        /* manipulation and information */
        
        void scale ( Float s ); 

        bool isfinite() const;
        
        FloatVector getDiagonal() const;

        MatrixCSR getTranspose() const;

        
        /* access and information to internal data */
        
        const int*   getA() const;
        
        const int*   getC() const;
        
        const Float* getV() const;

        int getnumberofentries() const;

        int getnumberofzeroentries() const;

        int getmaxrowwidth() const;

        Float eigenvalueupperbound() const;

        void compressentries() const;

        /* Memory size */
        
        std::size_t memorysize() const;
        
        
        
    private:

        std::vector<int>   A;
        std::vector<int>   C; // column index of each term 
        std::vector<Float> V; // numerical value of each term
    
};




MatrixCSR MatrixCSRMultiplication( const MatrixCSR& mat1, const MatrixCSR& mat2 );

MatrixCSR MatrixCSRMultiplication_reduced( const MatrixCSR& mat1, const MatrixCSR& mat2 );

MatrixCSR MatrixCSRAddition( const MatrixCSR& mat1, const MatrixCSR& mat2, Float s1, Float s2 );

DiagonalOperator InverseDiagonalPreconditioner( const MatrixCSR& mat );



inline MatrixCSR operator+( const MatrixCSR& mat1, const MatrixCSR& mat2 )
{
//     LOG << "ADD" << nl; 
    return MatrixCSRAddition( mat1, mat2, 1., 1. );
}

inline MatrixCSR operator-( const MatrixCSR& mat1, const MatrixCSR& mat2 )
{
    return MatrixCSRAddition( mat1, mat2, 1., -1. );
}

inline MatrixCSR operator&( const MatrixCSR& mat1, const MatrixCSR& mat2 )
{
    return MatrixCSRMultiplication_reduced( mat1, mat2 );
}

inline MatrixCSR operator*( Float s, const MatrixCSR& mat )
{
    auto foo = mat;
    foo.scale(s);
    return foo;
}

inline MatrixCSR operator*( const MatrixCSR& mat, Float s )
{
    return s * mat;
}




#endif
