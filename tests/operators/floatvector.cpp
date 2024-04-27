
#include "../../basic.hpp"
#include "../../operators/floatvector.hpp"


int main( int argc, char *argv[] )
{
        LOG << "Unit Test for Vector class" << nl;
        
        if(true)
        {
        
            FloatVector a(5);
            
            a.check();
            
            LOG << "Should be zero vector:" << nl;
            a.zero();
            LOG << a << nl;
            
            for( int i = 0; i < 5; i++ )
                    a.setentry( i, i+1 );
            LOG << "Should be ascending numbers:" << nl;
            LOG << a.data_as_text() << nl;
            
            LOG << "Should be the middle entries:" << nl;
            LOG << a.getslice(1,3).data_as_text() << nl;
            
            FloatVector b(a);
            LOG << "Should be the same again:" << nl;
            LOG << b.data_as_text() << nl;
            
            LOG << "Should be multiples of PI:" << nl;
            LOG << (3.141 * a).data_as_text() << nl;
            
            LOG << "Should be negative of original vector:" << nl;
            LOG << (-a).data_as_text() << nl;
            
            FloatVector t(5);
            
            for( int i = 0; i < 5; i++ )
                    t.setentry( i, 3. * i+1 );
            LOG << "Should be other ascending numbers: 3*( i + 1) " << nl;
            LOG << t.data_as_text() << nl;
            
            LOG << "Next the sum of two vectors:" << nl;
            LOG << (a+t).data_as_text() << nl;
            LOG << "Then their difference:" << nl;
            LOG << (a-t).data_as_text() << nl;
            
            LOG << "Now the scalar product with itself:" << nl;
            LOG << a*a << nl;
            
            LOG << "Copy the middle slice:" << nl;
            a.setslice( 1, t.getslice(1,3) );
            LOG << a.data_as_text() << nl;
            
            LOG << "Add the middle slice:" << nl;
            a.addslice( 1, t.getslice(1,3), 1000. );
            LOG << a.data_as_text() << nl;
            
            FloatVector e(0);
            LOG << "Should be the zero-dimensional vector:" << nl;
            LOG << e.data_as_text() << nl;
            
        }
        
        if(true)
        {
        
            FloatVector a { 1, 3, 0 };
            FloatVector b { -1, -3, 0 };
            FloatVector c { 2, 4, 1 };
            FloatVector d { -5, -4, -3 };
            FloatVector e { 0,0,0 };
            
            LOG << FloatVector {3} << nl;
            
            LOG << "positive:     (no ) " << a.ispositive() << nl;
            LOG << "negative:     (no ) " << a.isnegative() << nl;
            LOG << "non-negative: (yes) " << a.isnonnegative() << nl;
            LOG << "non-positive: (no ) " << a.isnonpositive() << nl;
            LOG << "zero:         (no ) " << a.iszero() << nl;
            
            LOG << "positive:     (no ) " << b.ispositive() << nl;
            LOG << "negative:     (no ) " << b.isnegative() << nl;
            LOG << "non-negative: (no ) " << b.isnonnegative() << nl;
            LOG << "non-positive: (yes) " << b.isnonpositive() << nl;
            LOG << "zero:         (no ) " << b.iszero() << nl;
            
            LOG << "positive:     (yes) " << c.ispositive() << nl;
            LOG << "negative:     (no ) " << c.isnegative() << nl;
            LOG << "non-positive: (no ) " << c.isnonpositive() << nl;
            LOG << "non-negative: (yes) " << c.isnonnegative() << nl;
            LOG << "zero:         (no ) " << c.iszero() << nl;
            
            LOG << "positive:     (no ) " << d.ispositive()   << nl;
            LOG << "negative:     (yes) " << d.isnegative() << nl;
            LOG << "non-positive: (yes) " << d.isnonpositive() << nl;
            LOG << "non-negative: (no ) " << d.isnonnegative() << nl;
            LOG << "zero:         (no ) " << d.iszero() << nl;
            
            LOG << "positive:     (no ) " << e.ispositive()   << nl;
            LOG << "negative:     (no ) " << e.isnegative()   << nl;
            LOG << "non-positive: (yes) " << e.isnonpositive() << nl;
            LOG << "non-negative: (yes) " << e.isnonnegative() << nl;
            LOG << "zero:         (yes) " << e.iszero()       << nl;
            
            LOG << "norm: (3.162)" << a.norm() << space << a.normalize().norm() << nl;
            LOG << "norm: (3.162)" << b.norm() << space << b.normalize().norm() << nl;
            LOG << "norm: (4.582)" << c.norm() << space << c.normalize().norm() << nl;
            LOG << "norm: (7.071)" << d.norm() << space << d.normalize().norm() << nl;
            LOG << "norm: (0.000)" << e.norm() << nl;
            
        }
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;

        return 0;
}
