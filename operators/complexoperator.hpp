#ifndef INCLUDEGUARD_OPERATOR_COMPLEXOPERATOR_HPP
#define INCLUDEGUARD_OPERATOR_COMPLEXOPERATOR_HPP

#include <utility>

template<typename T>
inline T* allocate_(T& t)
{
    return std::move(t).pointer_to_heir();
}

template<typename T>
inline T* allocate_(const T& )
{
    return nullptr;
}

#include "../basic.hpp"
#include "floatvector.hpp"
#include "linearoperator.hpp"
#include "simpleoperators.hpp"


/************************
****
****  Class for Complex Operators 
****  
************************/



FloatVector RealPart( const FloatVector& vec );
FloatVector ImaginaryPart( const FloatVector& vec );
FloatVector ComplexFloatVector( const FloatVector& real, const FloatVector& imag );

class ComplexOperator final
: public LinearOperator 
{

    public:

        /* Constructors */
        
        template<typename OpReal>
        explicit ComplexOperator( OpReal&& Real )
        : ComplexOperator( std::forward<OpReal>(Real), ZeroOperator(Real.getdimout(),Real.getdimin()) )
        {}
        
        template<typename OpReal, typename OpImag>
        ComplexOperator( OpReal&& Real, OpImag&& Imag )
        : LinearOperator( Real.getdimout(), Real.getdimin() ), 
        part_real( std::is_lvalue_reference<OpReal>::value ? &Real : allocate_<typename std::remove_reference<OpReal>::type>( Real ) ),
        part_imag( std::is_lvalue_reference<OpImag>::value ? &Imag : allocate_<typename std::remove_reference<OpImag>::type>( Imag ) ),
        managing_real( not std::is_lvalue_reference<OpReal>::value ),
        managing_imag( not std::is_lvalue_reference<OpImag>::value )
        {
            ComplexOperator::check();
        }
        
        ComplexOperator( int dim, Float real, Float imag )
        : ComplexOperator( ScalingOperator(dim,real), ScalingOperator(dim,imag) )
        {}
        

        /* standard methods for operators */
        
        ComplexOperator()                                        = delete;
        ComplexOperator( const ComplexOperator& )                = delete;
        ComplexOperator& operator=( const ComplexOperator& vec ) = delete;

        ComplexOperator( ComplexOperator&& op )
        : LinearOperator( op.getdimout(), op.getdimin() ), 
        part_real( op.part_real ), part_imag( op.part_imag ),
        managing_real( op.managing_real ), managing_imag( op.managing_imag )
        {
            op.part_real = nullptr; op.managing_real = false;
            op.part_imag = nullptr; op.managing_imag = false;
        }

        ComplexOperator& operator=( ComplexOperator&& op ) = delete;
        // {
        //     part_real = op.part_real; managing_part_real = op.managing_part_real; op.part_real = nullptr; op.managing_part_real = false;
        //     part_imag = op.part_imag; managing_part_imag = op.managing_part_imag; op.part_imag = nullptr; op.managing_part_imag = false;
        //     return *this;
        // }


        virtual ~ComplexOperator();
        
        /* standard interface */
        
        virtual void check() const override;
        
        virtual std::string text() const override;
        
        virtual std::string text( const bool embellish ) const;

        /* OTHER METHODS */
        
        virtual ComplexOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override;
        
    private:

        const LinearOperator* part_real;
        const LinearOperator* part_imag;
        bool managing_real;
        bool managing_imag;
            
};



  
  
  
  

#endif
