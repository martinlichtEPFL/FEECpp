
#include <utility>

#include "systemsolver.hpp"

#include "../operators/floatvector.hpp"


static const bool cppsys_restart_on_full_dimension = false;
static const bool cppsys_restart_before_finish     = false;



  
int BlockHerzogSoodhalterMethod( 
    FloatVector& x_A, 
    FloatVector& x_C, 
    const FloatVector& b_A, 
    const FloatVector& b_C, 
    const LinearOperator& A, const LinearOperator& Bt, const LinearOperator& B, const LinearOperator& C, 
    Float tolerance,
    int print_modulo,
    const LinearOperator& PAinv, const LinearOperator& PCinv
) {

    x_A.check();
    x_C.check();
    
    b_A.check();
    b_C.check();
    
    A.check(); C.check(); B.check(); Bt.check(); 

    assert( A.getdimin() == A.getdimout() );
    assert( C.getdimin() == C.getdimout() );
    
    const int dimension_A = A.getdimin();
    const int dimension_C = C.getdimin();

    assert( B.getdimin()  == dimension_A );
    assert( B.getdimout() == dimension_C );
    
    assert( Bt.getdimin()  == B.getdimout() );
    assert( Bt.getdimout() == B.getdimin()  );
    
    assert( b_A.getdimension() == dimension_A );
    assert( b_C.getdimension() == dimension_C );
    assert( x_A.getdimension() == dimension_A );
    assert( x_C.getdimension() == dimension_C );
    
    assert( PAinv.getdimin() == PAinv.getdimout() );
    assert( PCinv.getdimin() == PCinv.getdimout() );
    assert( dimension_A == PAinv.getdimout() );
    assert( dimension_C == PCinv.getdimout() );

    FloatVector v0_A( dimension_A, 0. );
    FloatVector v1_A( dimension_A, 0. );
    FloatVector w0_A( dimension_A, 0. );
    FloatVector w1_A( dimension_A, 0. );
    FloatVector  z_A( dimension_A, 0. );
    
    FloatVector v0_C( dimension_C, 0. );
    FloatVector v1_C( dimension_C, 0. );
    FloatVector w0_C( dimension_C, 0. );
    FloatVector w1_C( dimension_C, 0. );
    FloatVector  z_C( dimension_C, 0. );
    
    FloatVector vn_A( dimension_A, 0. );
    FloatVector wn_A( dimension_A, 0. );
    FloatVector zn_A( dimension_A, 0. );
    
    FloatVector vn_C( dimension_C, 0. );
    FloatVector wn_C( dimension_C, 0. );
    FloatVector zn_C( dimension_C, 0. );
    
    FloatVector  m_A( dimension_A, 0. );
    FloatVector  m_C( dimension_C, 0. );
    
    FloatVector p_A( dimension_A, 0. );
    FloatVector p_C( dimension_C, 0. );
    
    
    Float mu_A = notanumber;
    Float mu_C = notanumber;
    
    Float gamma = notanumber;
    
    Float eta   = notanumber;
    Float eta_A = notanumber;
    Float eta_C = notanumber;
    
    Float s0 = notanumber;
    Float s1 = notanumber;
    Float c0 = notanumber;
    Float c1 = notanumber;
    
    int max_iteration_count = dimension_A + dimension_C;
    int recent_iteration_count = 0;

    if( print_modulo >= 0 ) LOGPRINTF( "START Block Herzog-Soodhalter CSR\n" );

    while( recent_iteration_count < max_iteration_count ){
        
        bool restart_condition = ( recent_iteration_count == 0 ) or ( cppsys_restart_on_full_dimension and recent_iteration_count );;
        
        bool residual_seems_small = ( recent_iteration_count != 0 ) and ( absolute(eta) < tolerance );
        
        if( restart_condition or ( residual_seems_small and cppsys_restart_before_finish ) ) UNLIKELY {
            
            v0_A.zero(); w0_A.zero(); w1_A.zero();
            v0_C.zero(); w0_C.zero(); w1_C.zero();

            v1_A = b_A - A * x_A - Bt * x_C;
            z_A = PAinv * v1_A;
            v1_C = b_C - B * x_A - C  * x_C;
            z_C = PCinv * v1_C;
            
            gamma = sqrt( v1_A * z_A + v1_C * z_C );
            
            v1_A /= gamma;
            v1_C /= gamma;
            z_A /= gamma;
            z_C /= gamma;
            
            Float psi_A = z_A * v1_A; 
            Float psi_C = z_C * v1_C;
            mu_A = psi_A; 
            mu_C = psi_C; 
            
            
            m_A = v1_A;
            m_C = v1_C;

            eta = gamma; 
            eta_A = gamma * sqrt( psi_A );
            eta_C = gamma * sqrt( psi_C );
            
            assert( gamma >= 0. );
            
            s0 = s1 = 0;
            c0 = c1 = 1;
            
            if( print_modulo >= 0 ) {
                LOGPRINTF( "(%d/%d) RESTARTED: Residual norm is %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) absolute(eta), (long double)tolerance );
                LOGPRINTF( "(%d/%d)            Gamma: %.9Le Eta_A %.9Le Eta_C %.9Le\n", recent_iteration_count, max_iteration_count, (long double)gamma, (long double)eta_A, (long double)eta_C );
            }

        }
        
        bool residual_is_small = ( absolute(eta) < tolerance );
        
        if( residual_is_small ) UNLIKELY 
            break;

            
        {
            

            p_A = A * z_A + Bt * z_C;
            p_C = B * z_A + C  * z_C;
 
            Float delta = z_A * p_A + z_C * p_C;

            vn_A = p_A - delta * v1_A - gamma * v0_A;
            vn_C = p_C - delta * v1_C - gamma * v0_C;
            
            zn_A = PAinv * vn_A;
            zn_C = PCinv * vn_C;

            Float gamma_n = std::sqrt( zn_A * vn_A + zn_C * vn_C );
 
            vn_A /= gamma_n; zn_A /= gamma_n;
            vn_C /= gamma_n; zn_C /= gamma_n;

            Float alpha_0 = c1 * delta - c0 * s1 * gamma;
            assert( alpha_0 * alpha_0 + gamma_n * gamma_n > 0. );
            Float alpha_1 = std::sqrt( alpha_0 * alpha_0 + gamma_n * gamma_n );
            Float alpha_2 = s1 * delta + c0 * c1 * gamma;
            Float alpha_3 = s0 * gamma;
 
            assert( alpha_1 > 0. );

            Float cn = alpha_0 / alpha_1;
            Float sn = gamma_n / alpha_1;
            
            
            Float theta_A = m_A * zn_A;
            Float theta_C = m_C * zn_C; 
            
            Float psi_A = zn_A * vn_A;
            Float psi_C = zn_C * vn_C; 

            m_A = - sn * m_A + cn * vn_A;
            m_C = - sn * m_C + cn * vn_C;

            wn_A = ( z_A - alpha_3 * w0_A - alpha_2 * w1_A ) / alpha_1;
            wn_C = ( z_C - alpha_3 * w0_C - alpha_2 * w1_C ) / alpha_1;

            x_A = x_A + cn * eta * wn_A;
            x_C = x_C + cn * eta * wn_C;
 
            mu_A = sn * sn * mu_A - 2 * sn * cn * theta_A + cn * cn * psi_A;
            mu_C = sn * sn * mu_C - 2 * sn * cn * theta_C + cn * cn * psi_C;

            eta   = -sn * eta;
            eta_A = eta * sqrt( psi_A );
            eta_C = eta * sqrt( psi_C );

            v0_A = v1_A; v1_A = vn_A;
            v0_C = v1_C; v1_C = vn_C;
            w0_A = w1_A; w1_A = wn_A;
            w0_C = w1_C; w1_C = wn_C;
            
            z_A = zn_A;
            z_C = zn_C;
            
            
            gamma = gamma_n;
            
            c0 = c1; c1 = cn;
            s0 = s1; s1 = sn;

            FloatVector r_A = b_A - A * x_A - Bt * x_C;
            FloatVector r_C = b_C - B * x_A - C  * x_C;
            Float r = sqrt( r_A * r_A + r_C * r_C );

            bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
            
            if( print_modulo > 0 and print_condition ) {
                LOGPRINTF( "(%d/%d)   INTERIM: Full Residual norm is %.9Le\n", recent_iteration_count, max_iteration_count, (long double)r );
                LOGPRINTF( "(%d/%d)            eta_A=%.9Le eta_C=%.9Le\n", recent_iteration_count, max_iteration_count, (long double)eta_A, (long double)eta_C );
            }

        }

        
        Float recent_deviation = eta;
        
        /* Print information */
        
        bool print_condition = ( print_modulo > 0 and recent_iteration_count % print_modulo == 0 );
        
        if( print_modulo >= 0 and print_condition ) UNLIKELY {
            LOGPRINTF( "(%d/%d)   INTERIM: Residual norm is %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double) recent_deviation, (long double)tolerance );
            LOGPRINTF( "(%d/%d)            Gamma: %.9Le Eta: %.9Le\n", recent_iteration_count, max_iteration_count, (long double)gamma, (long double)eta );
        }
        
        recent_iteration_count++;
        
    }
    
    /* HOW DID WE FINISH ? */
    
    Float recent_deviation = absolute( eta );
        
    if( print_modulo >= 0 ) 
        LOGPRINTF( "(%d/%d)  FINISHED: Residual norm is %.9Le < %.9Le\n", recent_iteration_count, max_iteration_count, (long double)recent_deviation, (long double)tolerance );

    return recent_iteration_count;
}
  


