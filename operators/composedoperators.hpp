#ifndef INCLUDEGUARD_OPERATOR_COMPOSEDOPERATORS_HPP
#define INCLUDEGUARD_OPERATOR_COMPOSEDOPERATORS_HPP

#include <type_traits>
#include <utility>

    
template<typename T>
inline T* allocate(T& t)
{
    return std::move(t).pointer_to_heir();
}

template<typename T>
inline T* allocate(const T& )
{
    return nullptr;
}

#include "../basic.hpp"
#include "linearoperator.hpp"
#include "simpleoperators.hpp"


/************************
****
****  Classes for Composition of Operators 
****  - instantiates LinearOperator
****  
************************/







class RepeatedDiagonalBlockOperator final
: public LinearOperator 
{

    public:

        /* Constructors */
        /* standard methods for operators */
        /* standard interface */
        /* OTHER METHODS */
        
        RepeatedDiagonalBlockOperator()                                                     = delete;
        RepeatedDiagonalBlockOperator( const RepeatedDiagonalBlockOperator& )               = delete;
        RepeatedDiagonalBlockOperator& operator=( const RepeatedDiagonalBlockOperator& op ) = delete;
        
        RepeatedDiagonalBlockOperator( RepeatedDiagonalBlockOperator&& op )
        : LinearOperator( op.getdimout(), op.getdimin() ), 
        internal(op.internal), managing_internal(op.managing_internal), repetition( op.repetition )
        {
            op.internal = nullptr; op.managing_internal = false;
        }

        RepeatedDiagonalBlockOperator& operator=( RepeatedDiagonalBlockOperator&& op ) = delete;
        // {
        //     internal = op.internal; managing_internal = op.managing_internal; op.internal = nullptr; op.managing_internal = false;
        //     repetition = op.repetition;
        //     return *this;
        // }

        virtual RepeatedDiagonalBlockOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }

    public:
        
        template<typename Op>
        RepeatedDiagonalBlockOperator( Op&& op, int repetition )
        : LinearOperator( repetition * op.getdimout(), repetition * op.getdimin() ), 
        internal( std::is_lvalue_reference<Op>::value ? &op : allocate<typename std::remove_reference<Op>::type>( op ) ),
        managing_internal( not std::is_lvalue_reference<Op>::value ),
        repetition( repetition )
        {
            RepeatedDiagonalBlockOperator::check();
        }
        
        virtual ~RepeatedDiagonalBlockOperator() {
            if( internal != nullptr && managing_internal ) delete internal;
        }
        

        virtual void check() const override { 
            assert( repetition >= 0 );
            if( internal == nullptr ) return;
            internal->check();
            assert( getdimout() == repetition * internal->getdimout() );
            assert( getdimin()  == repetition * internal->getdimin()  );
        }
        
        virtual std::string text() const override
        {
            return "Repeated Diagonal Block Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + "\n"
                     + tab_each_line( internal->text() );
        }
        
        using LinearOperator::apply;
        virtual void apply( FloatVector& dest, const FloatVector& src, Float scaling ) const override {
            check();
            src.check();
            dest.check();
            
            assert( getdimin() == src.getdimension() );
            assert( getdimout() == dest.getdimension() );
            assert( src.getdimension()  % repetition == 0 );
            assert( dest.getdimension() % repetition == 0 );
            
            const int internal_dimout = internal->getdimout();
            const int internal_dimin  = internal->getdimin();
            
            for( int r = 0; r < repetition; r++ ) {
                auto src_slice = src.getslice( r * internal_dimin, internal_dimin );
                auto new_slice = internal->apply( src_slice, scaling );
                dest.setslice( r * internal_dimout, new_slice );
            }

        }

    private:

        const LinearOperator* internal;
        bool managing_internal;
        const int repetition;
    
};




































class ComposedOperator
: public LinearOperator 
{

    public:
    
        /* Constructors */
        /* standard methods for operators */
        /* standard interface */
        /* OTHER METHODS */
        
        ComposedOperator()                                        = delete;
        ComposedOperator( const ComposedOperator& )               = delete;
        ComposedOperator& operator=( const ComposedOperator& op ) = delete;
        
        ComposedOperator( ComposedOperator&& op )
        : LinearOperator( op.getdimout(), op.getdimin() ), 
        left(op.left), right(op.right), managing_left(op.managing_left), managing_right(op.managing_right)
        {
            op.left  = nullptr; op.managing_left  = false;
            op.right = nullptr; op.managing_right = false;
        }

        ComposedOperator& operator=( ComposedOperator&& op ) = delete;
        // {
        //     left  = op.left;  managing_left  = op.managing_left;  op.left = nullptr;  op.managing_left  = false;
        //     right = op.right; managing_right = op.managing_right; op.right = nullptr; op.managing_right = false;
        //     return *this;
        // }

    public:

        template<typename OpL, typename OpR>
        ComposedOperator( int dimout, int dimin, OpL&& L, OpR&& R )
        : LinearOperator( dimout, dimin ), 
        left(  std::is_lvalue_reference<OpL>::value ? &L : allocate<typename std::remove_reference<OpL>::type>( L ) ),
        right( std::is_lvalue_reference<OpR>::value ? &R : allocate<typename std::remove_reference<OpR>::type>( R ) ),
        managing_left(  not std::is_lvalue_reference<OpL>::value ),
        managing_right( not std::is_lvalue_reference<OpR>::value )
        {
            ComposedOperator::check();
        }
        
        virtual ~ComposedOperator() { 
            if( left  != nullptr && managing_left  ) delete left;
            if( right != nullptr && managing_right ) delete right;
        }

        virtual void check() const override { 
            LinearOperator::check(); 
            if( left != nullptr  ) left->check();
            if( right != nullptr ) right->check(); 
        }
        
        virtual std::string text() const override
        {
            return "Composed Operator" + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + "\n"
                     + tab_each_line( left->text() + '\n' + right->text() );
        }
        
    protected:

        const LinearOperator* left;
        const LinearOperator* right;
        bool managing_left;
        bool managing_right;
        
};








class ProduktOperator final
: public ComposedOperator 
{

    public:
    
        /* Constructors */
        /* standard methods for operators */
        /* standard interface */
        /* OTHER METHODS */
        
        ProduktOperator()                                       = delete;
        ProduktOperator( const ProduktOperator& )               = delete;
        ProduktOperator& operator=( const ProduktOperator& op ) = delete;
        
        ProduktOperator( ProduktOperator&& )                    = default;
        ProduktOperator& operator=( ProduktOperator&& op )      = default; 

        virtual ProduktOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }

    public:

        template<typename OpL, typename OpR>
        ProduktOperator( OpL&& L, OpR&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::forward<OpL>(L), std::forward<OpR>(R) )
        {
            ComposedOperator::check();
        }

        virtual ~ProduktOperator()
        {
//             ProduktOperator::check();
        }

        virtual void check() const override { 
            ComposedOperator::check();
            if( left == nullptr || right == nullptr ) assert( left == nullptr && right == nullptr );
            if( left == nullptr || right == nullptr ) return;
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimin() == right->getdimout() );
        }
        
        virtual std::string text() const override
        {
            return "Produkt Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + "\n"
                     + tab_each_line( left->text() + '\n' + right->text() );
        }
        
        using LinearOperator::apply;
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override {
            dest = scaling * left->apply( right->apply(add) );
        }

};

inline ProduktOperator operator*( const LinearOperator& left, const LinearOperator& right )
{ return ProduktOperator( left, right ); }

inline ProduktOperator operator*( const LinearOperator& left, LinearOperator&& right )
{ return ProduktOperator( left, std::move(right) ); }

inline ProduktOperator operator*( LinearOperator&& left, const LinearOperator& right )
{ return ProduktOperator( std::move(left), right ); }

inline ProduktOperator operator*( LinearOperator&& left, LinearOperator&& right )
{ return ProduktOperator( std::move(left), std::move(right) ); }

inline ProduktOperator operator-( LinearOperator&& op )
{
    return ProduktOperator( ScalingOperator( op.getdimout(), -1. ), std::move(op) );
}

inline ProduktOperator operator-( const LinearOperator& op )
{
    return ProduktOperator( ScalingOperator( op.getdimout(), -1. ), op );
}

inline ProduktOperator operator+( LinearOperator&& op )
{
    return ProduktOperator( IdentityOperator( op.getdimout() ), std::move(op) );
}

inline ProduktOperator operator+( const LinearOperator& op )
{
    return ProduktOperator( IdentityOperator( op.getdimout() ), op );
}

inline ProduktOperator operator*( Float s, LinearOperator&& op )
{
    return ProduktOperator( ScalingOperator( op.getdimout(), s ), std::move(op) );
}

inline ProduktOperator operator*( Float s, const LinearOperator& op )
{
    return ProduktOperator( ScalingOperator( op.getdimout(), s ), op );
}










class SummOperator final
: public ComposedOperator 
{

    public:
    
        /* Constructors */
        /* standard methods for operators */
        /* standard interface */
        /* OTHER METHODS */
        
        SummOperator()                                    = delete;
        SummOperator( const SummOperator& )               = delete;
        SummOperator& operator=( const SummOperator& op ) = delete;
        
        SummOperator( SummOperator&& )                    = default;
        SummOperator& operator=( SummOperator&& op )      = default;  

        virtual SummOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        }
        
    public:

        template<typename OpL, typename OpR>
        SummOperator( OpL&& L, OpR&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::forward<OpL>(L), std::forward<OpR>(R) )
        {
            ComposedOperator::check();
        }

        virtual ~SummOperator()
        {
//             SummOperator::check();
        }
        
        virtual void check() const override
        { 
            ComposedOperator::check();
            if( left == nullptr || right == nullptr ) assert( left == nullptr && right == nullptr );
            if( left == nullptr || right == nullptr ) return;
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == left->getdimin()  );
            assert( getdimout() == right->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimout() == right->getdimout() );
            assert( left->getdimin()  == right->getdimin()  );
        }
        
        virtual std::string text() const override
        {
            return "Summ Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + "\n"
                     + tab_each_line( left->text() + '\n' + right->text() );
        }
        
        
        using LinearOperator::apply;
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override
        {
            dest = scaling * left->apply( add ) + scaling * right->apply( add );
        }

};

inline SummOperator operator+( const LinearOperator& left, const LinearOperator& right )
{ return SummOperator( left, right ); }

inline SummOperator operator+( const LinearOperator& left, LinearOperator&& right )
{ return SummOperator( left, std::move(right) ); }

inline SummOperator operator+( LinearOperator&& left, const LinearOperator& right )
{ return SummOperator( std::move(left), right ); }

inline SummOperator operator+( LinearOperator&& left, LinearOperator&& right )
{ return SummOperator( std::move(left), std::move(right) ); }
 






class DiffOperator final
: public ComposedOperator 
{

    public:
    
        /* Constructors */
        /* standard methods for operators */
        /* standard interface */
        /* OTHER METHODS */
        
        DiffOperator()                                    = delete;
        DiffOperator( const DiffOperator& )               = delete;
        DiffOperator& operator=( const DiffOperator& op ) = delete;
        
        DiffOperator( DiffOperator&& )                    = default;
        DiffOperator& operator=( DiffOperator&& op )      = default; 

        virtual DiffOperator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        } 
        
    public:

        template<typename OpL, typename OpR>
        DiffOperator( OpL&& L, OpR&& R )
        : ComposedOperator( L.getdimout(), R.getdimin(), std::forward<OpL>(L), std::forward<OpR>(R) )
        {
            ComposedOperator::check();
        }
        
        virtual ~DiffOperator()
        {
//             DiffOperator::check();
        }
        
        virtual void check() const override
        { 
            ComposedOperator::check();
            if( left == nullptr || right == nullptr ) assert( left == nullptr && right == nullptr );
            if( left == nullptr || right == nullptr ) return;
            assert( getdimout() == left->getdimout() );
            assert( getdimin()  == left->getdimin()  );
            assert( getdimout() == right->getdimout() );
            assert( getdimin()  == right->getdimin()  );
            assert( left->getdimout() == right->getdimout() );
            assert( left->getdimin()  == right->getdimin()  );
        }
        
        virtual std::string text() const override
        {
            return "Diff Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + "\n"
                     + tab_each_line( left->text() + '\n' + right->text() );
        }
        
        
        using LinearOperator::apply;
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override
        {
            dest = scaling * left->apply( add ) - scaling * right->apply( add );
        }

};

inline DiffOperator operator-( const LinearOperator& left, const LinearOperator& right )
{ return DiffOperator( left, right ); }

inline DiffOperator operator-( const LinearOperator& left, LinearOperator&& right )
{ return DiffOperator( left, std::move(right) ); }

inline DiffOperator operator-( LinearOperator&& left, const LinearOperator& right )
{ return DiffOperator( std::move(left), right ); }

inline DiffOperator operator-( LinearOperator&& left, LinearOperator&& right )
{ return DiffOperator( std::move(left), std::move(right) ); }
 

 




























class Block2x2Operator
: public LinearOperator 
{

    public:
    
        /* Constructors */
        /* standard methods for operators */
        /* standard interface */
        /* OTHER METHODS */
        
        Block2x2Operator()                                     = delete;
        Block2x2Operator( const Block2x2Operator& )            = delete;
        Block2x2Operator& operator=( const Block2x2Operator& ) = delete;
        
        Block2x2Operator( Block2x2Operator&& op )
        : LinearOperator( op.getdimout(), op.getdimin() ), 
        upperleft(op.upperleft), upperright(op.upperright), 
        lowerleft(op.lowerleft), lowerright(op.lowerright), 
        managing_upperleft(op.managing_upperleft), managing_upperright(op.managing_upperright),
        managing_lowerleft(op.managing_lowerleft), managing_lowerright(op.managing_lowerright)
        {
            op.upperleft  = nullptr; op.managing_upperleft  = false;
            op.upperright = nullptr; op.managing_upperright = false;
            op.lowerleft  = nullptr; op.managing_lowerleft  = false;
            op.lowerright = nullptr; op.managing_lowerright = false;
        }

        Block2x2Operator& operator=( Block2x2Operator&& op ) = delete;
        // {
        //     upperleft  = op.upperleft;  managing_upperleft  = op.managing_upperleft;  op.upperleft = nullptr;  op.managing_upperleft  = false;
        //     upperright = op.upperright; managing_upperright = op.managing_upperright; op.upperright = nullptr; op.managing_upperright = false;
        //     lowerleft  = op.lowerleft;  managing_lowerleft  = op.managing_lowerleft;  op.lowerleft = nullptr;  op.managing_lowerleft  = false;
        //     lowerright = op.lowerright; managing_lowerright = op.managing_lowerright; op.lowerright = nullptr; op.managing_lowerright = false;
        //     return *this;
        // } 

        virtual Block2x2Operator* pointer_to_heir() && override
        {
            return new typename std::remove_reference<decltype(*this)>::type( std::move(*this) );
        } 

                
    public:

        
        template<typename OpUL, typename OpUR, typename OpLL, typename OpLR>
        Block2x2Operator( int dimout, int dimin, 
                          OpUL&&  UL, OpUR&&  UR, 
                          OpLL&&  LL, OpLR&&  LR 
                        )
        : LinearOperator( UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin() ), 
        upperleft(  std::is_lvalue_reference<OpUL>::value ? &UL : allocate<typename std::remove_reference<OpUL>::type>( UL ) ),
        upperright( std::is_lvalue_reference<OpUR>::value ? &UR : allocate<typename std::remove_reference<OpUR>::type>( UR ) ),
        lowerleft(  std::is_lvalue_reference<OpLL>::value ? &LL : allocate<typename std::remove_reference<OpLL>::type>( LL ) ),
        lowerright( std::is_lvalue_reference<OpLR>::value ? &LR : allocate<typename std::remove_reference<OpLR>::type>( LR ) ),
        managing_upperleft(  not std::is_lvalue_reference<OpUL>::value ),
        managing_upperright( not std::is_lvalue_reference<OpUR>::value ),
        managing_lowerleft(  not std::is_lvalue_reference<OpLL>::value ),
        managing_lowerright( not std::is_lvalue_reference<OpLR>::value )
        {
            assert( dimout >= 0 and dimin >= 0 );
            assert( dimout == UL.getdimout() + LL.getdimout() );
            assert( dimin  == UL.getdimin()  + UR.getdimin()  );
            assert( UL.getdimin()  == LL.getdimin()  );
            assert( UR.getdimin()  == LR.getdimin()  );
            assert( UL.getdimout() == UR.getdimout() );
            assert( LL.getdimout() == LR.getdimout() );
            Block2x2Operator::check();
        }
        
        template<typename OpUL, typename OpUR, typename OpLL, typename OpLR>
        Block2x2Operator( OpUL&&  UL, OpUR&&  UR, OpLL&&  LL, OpLR&&  LR )
        : Block2x2Operator( 
            UL.getdimout() + LL.getdimout(), UL.getdimin() + UR.getdimin(), 
            std::forward<OpUL>(UL), std::forward<OpUR>(UR), std::forward<OpLL>(LL), std::forward<OpLR>(LR) 
        ){}
        

        using LinearOperator::apply;
        void apply( FloatVector& dest, const FloatVector& add, Float scaling ) const override {
            
            assert( dest.getdimension() == upperleft->getdimout()  + lowerleft->getdimout()  );
            assert( dest.getdimension() == upperright->getdimout() + lowerright->getdimout() );
            assert(  add.getdimension() == upperleft->getdimin()   + upperright->getdimin()  );
            assert(  add.getdimension() == lowerleft->getdimin()   + lowerright->getdimin()  );
            
            
            const auto left  = add.getslice( 0,                     upperleft ->getdimin() );
            const auto right = add.getslice( upperleft->getdimin(), upperright->getdimin() );
            
            const auto upper = (*upperleft) * left + (*upperright) * right;
            const auto lower = (*lowerleft) * left + (*lowerright) * right;
            
            dest.setslice( 0,                    upper );
            dest.setslice( upper.getdimension(), lower );
            
            dest.scale( scaling );
        }

        
        
        
        virtual ~Block2x2Operator() {
            if( upperleft  != nullptr && managing_upperleft  ) delete upperleft;
            if( upperright != nullptr && managing_upperright ) delete upperright;
            if( lowerleft  != nullptr && managing_lowerleft  ) delete lowerleft;
            if( lowerright != nullptr && managing_lowerright ) delete lowerright;
        }

        virtual void check() const override { 
            LinearOperator::check(); 
            if( upperleft == nullptr || upperright == nullptr || lowerleft == nullptr || lowerright == nullptr )
                assert( upperleft == nullptr && upperright == nullptr && lowerleft == nullptr && lowerright == nullptr );
            if( upperleft == nullptr || upperright == nullptr || lowerleft == nullptr || lowerright == nullptr ) return;
            upperleft->check();
            upperright->check(); 
            lowerleft->check();
            lowerright->check(); 
        }
        
        virtual std::string text() const override
        {
            return "Block Operator " + std::to_string(getdimout()) + "x" + std::to_string(getdimin()) + "\n"
                    + 
                    tab_each_line( 
                        upperleft->text() + '\n' + upperright->text() + '\n' + lowerleft->text() + '\n' + lowerright->text()
                    );
        }
        
    private:

        const LinearOperator* upperleft;
        const LinearOperator* upperright;
        const LinearOperator* lowerleft;
        const LinearOperator* lowerright;
        bool managing_upperleft;
        bool managing_upperright;
        bool managing_lowerleft;
        bool managing_lowerright;
        
};



















#endif
