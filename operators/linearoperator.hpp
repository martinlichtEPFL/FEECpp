#ifndef INCLUDEGUARD_OPERATOR_LINEAROPERATOR_HPP
#define INCLUDEGUARD_OPERATOR_LINEAROPERATOR_HPP




//#include <memory>
// #include <ostream>

#include "../basic.hpp"
#include "floatvector.hpp"



/******************
*** 
***  Linear Operator on an abstract level 
***  - applyAdd is a purely virtual function that implements y <- A x with optional scaling of the result
***  - all other apply adds are defined in terms of that, but overwriting in derived classes is possible
***
******************/





class LinearOperator
{

    public:
        
        /* Constructors */
        
        explicit LinearOperator( int );
        explicit LinearOperator( int, int );

        /* standard methods for operators */
        
                 LinearOperator()                             = delete;
        explicit LinearOperator( const LinearOperator& )      = default;
        explicit LinearOperator( LinearOperator&& )           = default;
        LinearOperator& operator=( const LinearOperator& op ) = default;
        LinearOperator& operator=( LinearOperator&& op )      = default;
        virtual ~LinearOperator();

        /* standard interface */
        
        virtual void check() const;

        virtual std::string text() const = 0;
        
        void print( std::ostream& os ) const;

        // // void lg() const { LOG << text() << nl; };
        
        /* OTHER METHODS */
        
        virtual LinearOperator* pointer_to_heir() && = 0;        
        
        
        int getdimin() const;

        int getdimout() const;
        
        bool issquare() const;
        
        /* Apply the operator */
        
        /* x := s A y */
        virtual FloatVector apply( const FloatVector& src, Float scaling = 1. ) const;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling = 1. ) const = 0;
        
        FloatVector createinputvector() const;
        FloatVector createoutputvector() const;
        
        
    private:
        
        int dimout;
        int dimin;
    
};
  
  
  
inline FloatVector operator*( const LinearOperator& op, const FloatVector& vec )
{
    op.check();
    vec.check();
    FloatVector ret( op.getdimout() );
    op.apply( ret, vec );
    vec.check();
    return ret;
}
  
template<typename Stream>
inline Stream& operator<<( Stream&& os, const LinearOperator& op )
{
    op.check();
    os << op.text(); // op.print( os );
    op.check();
    return os;
}
  


  
  
  
  
#endif
