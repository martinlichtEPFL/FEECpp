#ifndef INCLUDEGUARD_SPARSE_SPARSEMATRIX_HPP
#define INCLUDEGUARD_SPARSE_SPARSEMATRIX_HPP

#include <utility>
#include <vector>

#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "../operators/simpleoperators.hpp"

class DenseMatrix;


/************************
****
****  Class for Sparse Matrices  
****  - instantiates LinearOperator
****  
************************/




class SparseMatrix:
public LinearOperator /* every matrix is a linear operator */
{

    public:

        struct PACKED MatrixEntry
        {
            MatrixEntry(int row, int column, Float value) : row(row), column(column), value(value) {}
            MatrixEntry(){}
            int row;
            int column;
            Float value;
        };

        static_assert( sizeof(MatrixEntry) == 2 * sizeof(int) + sizeof(Float), "MatrixEntry takes too much memory" );

        enum class MatrixEntrySorting : unsigned char {
            rowwise,
            columnwise
        };

        /* Constructors */
        
        explicit SparseMatrix( int dimout, int dimin, int numentries = 0, 
                               std::function<MatrixEntry(int)> generator = [](int)->MatrixEntry{ return MatrixEntry(0,0,notanumber); } ); 
        // explicit SparseMatrix( int dimout, int dimin );
        explicit SparseMatrix( int dimout, int dimin, const std::vector<MatrixEntry>& entries );
        explicit SparseMatrix( int dimout, int dimin, const std::initializer_list<MatrixEntry>& entries );
        explicit SparseMatrix( const FloatVector& diagonal );
        explicit SparseMatrix( const ScalingOperator& matrix );
        explicit SparseMatrix( const DiagonalOperator& matrix );
        explicit SparseMatrix( const DenseMatrix& matrix );
        
        /* standard interface */ 
        
        SparseMatrix() = delete;
        SparseMatrix( const SparseMatrix& );
        SparseMatrix& operator=( const SparseMatrix& );
        SparseMatrix( SparseMatrix&& );
        SparseMatrix& operator=( SparseMatrix&& );
        virtual ~SparseMatrix();

        /* standard methods for operators */
        
        virtual void check() const override;
        virtual std::string text() const override;
        // virtual void printplain( std::ostream& ) const;


        /* OTHER METHODS */
        
        virtual SparseMatrix* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override;

        
        /* manipulation and information */
        
        void scale( Float s );

        bool isfinite() const;
        
        FloatVector getDiagonal() const;
        
        int getnumberofzeroentries() const;
        
        
        /* access and information to internal data */
        
        const std::vector<MatrixEntry>& getentries() const;
        
        std::vector<MatrixEntry>& getentries();
        
        int getnumberofentries() const;
        
        
        /* sorting entries */
        
        bool is_sorted( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;
        const SparseMatrix& sortentries( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;
        const SparseMatrix& sortandcompressentries( MatrixEntrySorting manner = MatrixEntrySorting::rowwise ) const;

        /* specific entry manipulations */
        
        void reserve( int ) const;
                
        const MatrixEntry& getentry( int ) const;
        
        MatrixEntry& getentry( int );
        
        void setentry( int, int, int, Float );
        
        void setentry( int, MatrixEntry );
        
        void appendentry( int, int, Float );
        
        void appendentry( MatrixEntry );
        
        void clearentries();
        
        /* obtain a transpose */

        SparseMatrix getTranspose() const;

        /* Memory size */
        
        std::size_t memorysize() const;
        
    private:

        mutable std::vector<MatrixEntry> entries; 
    
};



SparseMatrix SparseMatrixMultiplication( const SparseMatrix& left, const SparseMatrix& right );

// FloatVector InverseDiagonalPreconditioner( const SparseMatrix& mat );

DiagonalOperator InverseDiagonalPreconditioner( const SparseMatrix& mat );



inline SparseMatrix Transpose( const SparseMatrix& op )
{
    return op.getTranspose();
}

inline SparseMatrix operator&( const SparseMatrix& left, const SparseMatrix& right )
{
    return SparseMatrixMultiplication( left, right );
}

inline SparseMatrix operator*( Float s, const SparseMatrix& mat )
{
    auto foo = mat;
    foo.scale(s);
    return foo;
}

inline SparseMatrix operator*( const SparseMatrix& mat, Float s )
{
    return s * mat;
}

Float norm_sq_of_vector( const SparseMatrix& A, const FloatVector& vec );










  

#endif
