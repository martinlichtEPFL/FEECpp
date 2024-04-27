
#include "../basic.hpp"
#include "../operators/floatvector.hpp"
#include "../operators/linearoperator.hpp"
#include "sparsesolver.hpp"
#include "iterativesolver.hpp"
#include "../utility/random.hpp"

#include "nullspace.hpp"

std::vector<FloatVector> computeNullspace(
    const LinearOperator& SystemMatrix,
    const LinearOperator& MassMatrix,
    const int max_number_of_candidates,
    //
    const Float tolerance_residual, 
    const Float tolerance_zero
) {

    std::vector<FloatVector> nullvectorgallery;
    
    for( int no_candidate = 0; no_candidate < max_number_of_candidates; no_candidate++ )
    {
        
        /* 1. create a new candidate */
        
        LOG << "New candidate" << nl;

        FloatVector candidate( SystemMatrix.getdimin(), 0. ); 
        candidate.random(); 
        candidate.normalize( MassMatrix );
        
        /* 2. reduce the candidate to its nullspace component */
        /*    and orthogonalize against previous nullvectors  */
        
        // filter out non-nullspace of the candidate 
        
        {
            FloatVector rhs( SystemMatrix.getdimin(), 0. );
            FloatVector res( SystemMatrix.getdimin(), 0. );
        
            ConjugateResidualMethod solver( SystemMatrix );
            solver.print_modulo        = -1; //candidate.getdimension();
            solver.max_iteration_count = candidate.getdimension();
            solver.tolerance           = tolerance_residual;

            solver.solve( candidate, rhs );
            
            res = rhs - SystemMatrix * candidate;
            
            LOG << "\t\t\t (eucl) res:       " << res.norm() << nl;
            LOG << "\t\t\t (mass) res:       " << res.norm( MassMatrix ) << nl;
            LOG << "\t\t\t (eucl) x:         " << candidate.norm() << nl;
            LOG << "\t\t\t (mass) x:         " << candidate.norm( MassMatrix ) << nl;
        }
        
        // gram schmidt cycles 

        for( const auto& nullvector : nullvectorgallery ) {
            Float alpha = ( MassMatrix * candidate * nullvector ) / ( MassMatrix * nullvector * nullvector );
            LOG << alpha << nl;
            candidate = candidate - alpha * nullvector;
        }
        
        LOG << "\t\t\t (eucl) res:       " << (SystemMatrix * candidate).norm() << nl;
        LOG << "\t\t\t (mass) res:       " << (SystemMatrix * candidate).norm( MassMatrix ) << nl;
        LOG << "\t\t\t (eucl) x:         " << candidate.norm() << nl;
        LOG << "\t\t\t (mass) x:         " << candidate.norm( MassMatrix ) << nl;
        
        /* 3. check whether anything is left at all */
        /*    exit if below tolerance               */

        Float purified_mass = candidate.norm( MassMatrix );
        LOG << "\t\t\t Filtered mass: " << purified_mass << nl;
        
        if( purified_mass < tolerance_zero )
            break;
        
        /* 4. If enough left, normalize and add to gallery */            
        
        candidate.normalize( MassMatrix );

        LOG << "\t\t\t (norm eucl) x:  " << candidate.norm() << nl;
        LOG << "\t\t\t (norm mass) x:  " << candidate.norm( MassMatrix ) << nl;
        LOG << "\t\t\t (norm eucl) Ax: " << ( SystemMatrix* candidate ).norm() << nl;
        LOG << "\t\t\t (norm mass) Ax: " << ( SystemMatrix * candidate ).norm( MassMatrix ) << nl;

        LOG << "Accept vector: " << nullvectorgallery.size() + 1 << nl;

        nullvectorgallery.push_back( candidate );
    }
    
    
    
    LOG << "How much nullspace are our vectors? (" << tolerance_residual << ")" << nl;
    for( const auto& nullvector : nullvectorgallery ) {
        LOGPRINTF( "% 10.5Le\t", (long double)( SystemMatrix * nullvector ).norm( MassMatrix ) );
    }
    LOG << nl;
    
    LOG << "How orthonormal are our vectors?" << nl;
    for( const auto& nullvector1 : nullvectorgallery ) {
        for( const auto& nullvector2 : nullvectorgallery ) {
            Float mass_prod = ( MassMatrix * nullvector1) * nullvector2;
            LOGPRINTF( "% 10.5Le\t", (long double)mass_prod );
        }
        LOG << nl;
    }

    return nullvectorgallery;

}

