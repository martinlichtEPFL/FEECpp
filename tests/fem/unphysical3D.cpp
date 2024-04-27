

/**/

#include "../../basic.hpp"
#include "../../mesh/mesh.simplicial3D.hpp"
#include "../../mesh/examples3D.hpp"
#include "../../fem/local.polynomialmassmatrix.hpp"
#include "../../fem/global.massmatrix.hpp"
#include "../../fem/global.elevation.hpp"
#include "../../fem/utilities.hpp"
#include "../../fem/global.unphysical.hpp"
#include "../../utility/convergencetable.hpp"


using namespace std;

int main( int argc, char *argv[] )
{
        
        LOG << "Unit Test: (3D) degree elevation of interpolation preserves mass" << nl;
        
        // LOG << std::setprecision(10);

        LOG << "Initial mesh..." << nl;
        
        MeshSimplicial3D M = UnitSimplex3D();
        
        M.check();
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_scalar_field;
        std::vector<Float>                                          experiments_scalar_value;
        

        
        
        experiments_scalar_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                assert( vec.getdimension() == 3 );
                return FloatVector({ std::exp( vec[0] + vec[1] + vec[2] ) });
            }
        );

        experiments_scalar_value.push_back( 3.19452804946532511361521373028750 * 3.19452804946532511361521373028750 * 3.19452804946532511361521373028750 );


        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_vector_field;
        std::vector<Float>                                          experiments_vector_value;
        

        
        
        experiments_vector_field.push_back( 
            [](const FloatVector& vec) -> FloatVector{
                    assert( vec.getdimension() == 3 );
                    Float x = vec[0]; Float y = vec[1]; Float z = vec[2]; 
                    return FloatVector({ std::sin( x*y ), std::cos( z ), std::exp(x+y+z) });
                }
            );

        experiments_vector_value.push_back( 33.42616007376754121867328484399462245509699700639766495311 );
        
        
        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_pseudo_field;
        std::vector<Float>                                          experiments_pseudo_value;

        experiments_pseudo_field = experiments_vector_field;
        experiments_pseudo_value = experiments_vector_value;
        
        
        
        
        
        
        
        
        std::vector<std::function<FloatVector(const FloatVector&)>> experiments_volume_field;
        std::vector<Float>                                          experiments_volume_value;

        experiments_volume_field = experiments_scalar_field;
        experiments_volume_value = experiments_scalar_value;

        
        
        const int r_min = 0;
        
        const int r_max = 2;
        
        const int l_min = 0;
        
        const int l_max = 2;
        
        Float errors_scalar[ experiments_scalar_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_vector[ experiments_vector_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_pseudo[ experiments_pseudo_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        Float errors_volume[ experiments_volume_field.size() ][ l_max - l_min + 1 ][ r_max - r_min + 1 ];
        
        
        
        for( int l = 0; l < l_min; l++ )
            M.uniformrefinement();

        for( int l = l_min; l <= l_max; l++ ){
            
            LOG << "Level:" << space << l_min << " <= " << l << " <= " << l_max << nl;
            
            for( int r      = r_min; r      <=      r_max; r++      ) 
            {
                
                LOG << "Polydegree:" << space << r_min << " <= " << r << " <= " << r_max << nl;
                
                LOG << "assemble mass matrices..." << nl;
                
                SparseMatrix massmatrix_scalar = FEECBrokenMassMatrix( M, M.getinnerdimension(), 0, r );
                SparseMatrix massmatrix_vector = FEECBrokenMassMatrix( M, M.getinnerdimension(), 1, r );
                SparseMatrix massmatrix_pseudo = FEECBrokenMassMatrix( M, M.getinnerdimension(), 2, r );
                SparseMatrix massmatrix_volume = FEECBrokenMassMatrix( M, M.getinnerdimension(), 3, r );

                SparseMatrix canonical_scalar = FEECCanonicalizeBroken( M, M.getinnerdimension(), 0, r );
                SparseMatrix canonical_vector = FEECCanonicalizeBroken( M, M.getinnerdimension(), 1, r );
                SparseMatrix canonical_pseudo = FEECCanonicalizeBroken( M, M.getinnerdimension(), 2, r );
                SparseMatrix canonical_volume = FEECCanonicalizeBroken( M, M.getinnerdimension(), 3, r );
                
                SparseMatrix randomize_scalar = FEECRandomizeBroken( M, M.getinnerdimension(), 0, r, notanumber );
                SparseMatrix randomize_vector = FEECRandomizeBroken( M, M.getinnerdimension(), 1, r, notanumber );
                SparseMatrix randomize_pseudo = FEECRandomizeBroken( M, M.getinnerdimension(), 2, r, notanumber );
                SparseMatrix randomize_volume = FEECRandomizeBroken( M, M.getinnerdimension(), 3, r, notanumber );
                
                // LOG << massmatrix_scalar << nl;
                // LOG << massmatrix_vector << nl;
                // LOG << massmatrix_pseudo << nl;
                // LOG << massmatrix_volume << nl;

                for( int i = 0; i < experiments_scalar_field.size(); i++ ){

                    const auto& scalarfield = experiments_scalar_field[i];
        
                    FloatVector interpol       = Interpolation( M, M.getinnerdimension(), 0, r, scalarfield );

                    FloatVector interpol_canon = canonical_scalar * interpol; // FEECCanonicalizeBroken( 3, 0, r, interpol );

                    FloatVector interpol_randm = randomize_scalar * interpol; // FEECRandomizeBroken( 3, 0, r, interpol );

                    Float mass_ii = absolute( interpol * ( massmatrix_scalar * interpol ) );

                    Float mass_ci = absolute( interpol_canon * ( massmatrix_scalar * interpol ) );

                    Float mass_cc = absolute( interpol_canon * ( massmatrix_scalar * interpol_canon ) );

                    Float mass_rr = absolute( interpol_randm * ( massmatrix_scalar * interpol_randm ) );

                    Float mass_cr = absolute( interpol_canon * ( massmatrix_scalar * interpol_randm ) );

                    assert( mass_ii >= -desired_closeness );
                    assert( mass_ci >= -desired_closeness );
                    assert( mass_cc >= -desired_closeness );
                    assert( mass_rr >= -desired_closeness );
                    assert( mass_cr >= -desired_closeness );
                    
                    Float error_mass = maxabs( mass_ii - mass_ci, mass_ii - mass_cc, mass_ii - mass_rr, mass_ii - mass_cr );
                    
                    LOG << mass_ii - mass_ci << space << mass_ii - mass_cc << space << mass_ii - mass_rr << space << mass_ii - mass_cr << nl;
                    
                    errors_scalar[i][l-l_min][r-r_min] = error_mass;
                    
                }
                
                for( int i = 0; i < experiments_vector_field.size(); i++ ){

                    const auto& vectorfield = experiments_vector_field[i];
        
                    FloatVector interpol       = Interpolation( M, M.getinnerdimension(), 1, r, vectorfield );

                    FloatVector interpol_canon = canonical_vector * interpol; // FEECCanonicalizeBroken( 3, 1, r, interpol );

                    FloatVector interpol_randm = randomize_vector * interpol; // FEECRandomizeBroken( 3, 1, r, interpol );

                    Float mass_ii = absolute( interpol * ( massmatrix_vector * interpol ) );

                    Float mass_ci = absolute( interpol_canon * ( massmatrix_vector * interpol ) );

                    Float mass_cc = absolute( interpol_canon * ( massmatrix_vector * interpol_canon ) );

                    Float mass_rr = absolute( interpol_randm * ( massmatrix_vector * interpol_randm ) );

                    Float mass_cr = absolute( interpol_canon * ( massmatrix_vector * interpol_randm ) );

                    assert( mass_ii >= -desired_closeness );
                    assert( mass_ci >= -desired_closeness );
                    assert( mass_cc >= -desired_closeness );
                    assert( mass_rr >= -desired_closeness );
                    assert( mass_cr >= -desired_closeness );
                    
                    Float error_mass = maxabs( mass_ii - mass_ci, mass_ii - mass_cc, mass_ii - mass_rr, mass_ii - mass_cr );
                    
                    LOG << mass_ii - mass_ci << space << mass_ii - mass_cc << space << mass_ii - mass_rr << space << mass_ii - mass_cr << nl;
                    
                    errors_vector[i][l-l_min][r-r_min] = error_mass;
                    
                }
                
                for( int i = 0; i < experiments_pseudo_field.size(); i++ ){

                    const auto& pseudofield = experiments_pseudo_field[i];
        
                    FloatVector interpol       = Interpolation( M, M.getinnerdimension(), 2, r, pseudofield );

                    FloatVector interpol_canon = canonical_pseudo * interpol; // FEECCanonicalizeBroken( 3, 2, r, interpol );

                    FloatVector interpol_randm = randomize_pseudo * interpol; // FEECRandomizeBroken( 3, 2, r, interpol );

                    Float mass_ii = absolute( interpol * ( massmatrix_pseudo * interpol ) );

                    Float mass_ci = absolute( interpol_canon * ( massmatrix_pseudo * interpol ) );

                    Float mass_cc = absolute( interpol_canon * ( massmatrix_pseudo * interpol_canon ) );

                    Float mass_rr = absolute( interpol_randm * ( massmatrix_pseudo * interpol_randm ) );

                    Float mass_cr = absolute( interpol_canon * ( massmatrix_pseudo * interpol_randm ) );

                    assert( mass_ii >= -desired_closeness );
                    assert( mass_ci >= -desired_closeness );
                    assert( mass_cc >= -desired_closeness );
                    assert( mass_rr >= -desired_closeness );
                    assert( mass_cr >= -desired_closeness );
                    
                    Float error_mass = maxabs( mass_ii - mass_ci, mass_ii - mass_cc, mass_ii - mass_rr, mass_ii - mass_cr );
                    
                    LOG << mass_ii - mass_ci << space << mass_ii - mass_cc << space << mass_ii - mass_rr << space << mass_ii - mass_cr << nl;
                    
                    errors_pseudo[i][l-l_min][r-r_min] = error_mass;
                    
                }
                
                for( int i = 0; i < experiments_volume_field.size(); i++ ){

                    const auto& volumefield = experiments_volume_field[i];
        
                    FloatVector interpol       = Interpolation( M, M.getinnerdimension(), 3, r, volumefield );

                    FloatVector interpol_canon = canonical_volume * interpol; // FEECCanonicalizeBroken( 3, 3, r, interpol );

                    FloatVector interpol_randm = randomize_volume * interpol; // FEECRandomizeBroken( 3, 3, r, interpol );

                    Float mass_ii = absolute( interpol * ( massmatrix_volume * interpol ) );

                    Float mass_ci = absolute( interpol_canon * ( massmatrix_volume * interpol ) );

                    Float mass_cc = absolute( interpol_canon * ( massmatrix_volume * interpol_canon ) );

                    Float mass_rr = absolute( interpol_randm * ( massmatrix_volume * interpol_randm ) );

                    Float mass_cr = absolute( interpol_canon * ( massmatrix_volume * interpol_randm ) );

                    assert( mass_ii >= -desired_closeness );
                    assert( mass_ci >= -desired_closeness );
                    assert( mass_cc >= -desired_closeness );
                    assert( mass_rr >= -desired_closeness );
                    assert( mass_cr >= -desired_closeness );
                    
                    Float error_mass = maxabs( mass_ii - mass_ci, mass_ii - mass_cc, mass_ii - mass_rr, mass_ii - mass_cr );
                    
                    LOG << mass_ii - mass_ci << space << mass_ii - mass_cc << space << mass_ii - mass_rr << space << mass_ii - mass_cr << nl;
                    
                    errors_volume[i][l-l_min][r-r_min] = error_mass;
                    
                }
                
            }
            
            if( l != l_max )
            {
                LOG << "Refinement..." << nl;
            
                M.uniformrefinement();

                M.shake_interior_vertices();
            }
            
        } 
    
        LOG << "Convergence tables" << nl;
    
        ConvergenceTable contable_scalar[ experiments_scalar_field.size() ];
        ConvergenceTable contable_vector[ experiments_vector_field.size() ];
        ConvergenceTable contable_pseudo[ experiments_vector_field.size() ];
        ConvergenceTable contable_volume[ experiments_volume_field.size() ];
        
        for( int r = r_min; r <= r_max; r++ ) 
        {
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i].table_name = "Rounding errors scalar E" + std::to_string(i);
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i].table_name = "Rounding errors vector E" + std::to_string(i);
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                contable_pseudo[i].table_name = "Rounding errors pseudo E" + std::to_string(i);
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i].table_name = "Rounding errors volume E" + std::to_string(i);

            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                contable_scalar[i] << printf_into_string("R%d", r-r_min );
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                contable_vector[i] << printf_into_string("R%d", r-r_min );
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                contable_pseudo[i] << printf_into_string("R%d", r-r_min );
            for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                contable_volume[i] << printf_into_string("R%d", r-r_min );        
        }
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
        for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i] << nl; 
        for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
        
        for( int l = l_min; l <= l_max; l++ ) 
        {
            
            for( int r = r_min; r <= r_max; r++ ) 
            {
                
                for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                    contable_scalar[i] << errors_scalar[i][l-l_min][r-r_min];
            
                for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                    contable_vector[i] << errors_vector[i][l-l_min][r-r_min];
            
                for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                    contable_pseudo[i] << errors_pseudo[i][l-l_min][r-r_min];
            
                for( int i = 0; i < experiments_volume_field.size(); i++ ) 
                    contable_volume[i] << errors_volume[i][l-l_min][r-r_min];
            
            }

            for( int i = 0; i < experiments_scalar_field.size(); i++ ) contable_scalar[i] << nl; 
            for( int i = 0; i < experiments_vector_field.size(); i++ ) contable_vector[i] << nl; 
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) contable_pseudo[i] << nl; 
            for( int i = 0; i < experiments_volume_field.size(); i++ ) contable_volume[i] << nl; 
            
        }

        
        
        
        
        LOG << "Convergence tables: scalars" << nl;
        for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
        {
            contable_scalar[i].lg(); 
            LOG << "-------------------" << nl;
        }
        
        LOG << "Convergence tables: vectors" << nl;
        for( int i = 0; i < experiments_vector_field.size(); i++ ) 
        {
            contable_vector[i].lg(); 
            LOG << "-------------------" << nl;
        }
        
        LOG << "Convergence tables: pseudos" << nl;
        for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
        {
            contable_pseudo[i].lg(); 
            LOG << "-------------------" << nl;
        }
        
        LOG << "Convergence tables: volumes" << nl;
        for( int i = 0; i < experiments_volume_field.size(); i++ )
        {
            contable_volume[i].lg(); 
            LOG << "-------------------" << nl;
        }
        
        
        
        
        
        LOG << "Check that differences are small: " << desired_closeness << nl;
        
        for( int l      = l_min; l      <=      l_max; l++      ) 
        for( int r      = r_min; r      <=      r_max; r++      ) 
        {
            for( int i = 0; i < experiments_scalar_field.size(); i++ ) 
                Assert( errors_scalar[i][l-l_min][r-r_min] < desired_closeness, errors_scalar[i][l-l_min][r-r_min], desired_closeness );
            
            for( int i = 0; i < experiments_vector_field.size(); i++ ) 
                Assert( errors_vector[i][l-l_min][r-r_min] < desired_closeness, errors_vector[i][l-l_min][r-r_min], desired_closeness );
            
            for( int i = 0; i < experiments_pseudo_field.size(); i++ ) 
                Assert( errors_pseudo[i][l-l_min][r-r_min] < desired_closeness, errors_pseudo[i][l-l_min][r-r_min], desired_closeness );
            
            for( int i = 0; i < experiments_volume_field.size(); i++ )
                Assert( errors_volume[i][l-l_min][r-r_min] < desired_closeness, errors_volume[i][l-l_min][r-r_min], desired_closeness );
        }
            
        
        LOG << "Finished Unit Test: " << ( argc > 0 ? argv[0] : "----" ) << nl;
        
        return 0;
}
