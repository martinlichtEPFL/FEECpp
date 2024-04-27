
#include <vector>
#include <set>

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "../dense/densematrix.hpp"
#include "../dense/functions.hpp"
#include "../sparse/sparsematrix.hpp"
#include "../sparse/matcsr.hpp"
#include "../mesh/mesh.hpp"

#include "../fem/lagrangematrices.hpp"
#include "../fem/utilities.hpp"
#include "../fem/polynomialmassmatrix.hpp"



SparseMatrix LagrangeBrokenMassMatrix( const Mesh& mesh, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_volumes = mesh.count_simplices( n );
    
    const int dim_in  = num_volumes * (n+1);
    const int dim_out = num_volumes * (n+1);
    
    
    // Set up sparse matrix
    
    SparseMatrix ret( dim_out, dim_in, num_volumes * (n+1)*(n+1) );
    
    
    // go over the n simplices and their k subsimplices
    for( int t  = 0; t  <  num_volumes; t++  )
    for( int v1 = 0; v1 <= n;           v1++ )
    for( int v2 = 0; v2 <= n;           v2++ )
    {
        
        SparseMatrix::MatrixEntry entry;
        
        int index_of_entry = t * (n+1)*(n+1) + v1 * (n+1) + v2;
        
        entry.row    = t * (n+1) + v1;
        
        entry.column = t * (n+1) + v2;
        
        Float measure = mesh.getMeasure( n, t );
        
        if( v1 == v2 )
            entry.value = 2. * factorial_numerical(n) * measure / factorial_numerical( 2 + n );
        else
            entry.value =      factorial_numerical(n) * measure / factorial_numerical( 2 + n );
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}




SparseMatrix LagrangeMassMatrix( const Mesh& mesh, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_volumes = mesh.count_simplices( n );
        
    const int num_vertices = mesh.count_simplices( 0 );
    
    const int dim_in  = num_vertices;
    const int dim_out = num_vertices;
    
    
    // Set up sparse matrix
    
    SparseMatrix ret( dim_out, dim_in, num_volumes * (n+1)*(n+1) );
    
    
    // go over the n simplices and their k subsimplices
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int t  = 0; t  <  num_volumes; t++  )
    {
        
        Float measure = mesh.getMeasure( n, t );
        
        for( int v1 = 0; v1 <= n;           v1++ )
        for( int v2 = 0; v2 <= n;           v2++ )
        {
            
            SparseMatrix::MatrixEntry entry;
            
            int index_of_entry = t * (n+1)*(n+1) + v1 * (n+1) + v2;
            
            int vertex1 = mesh.get_subsimplex( n, 0, t, v1 );
            int vertex2 = mesh.get_subsimplex( n, 0, t, v2 );
            
            entry.row    = mesh.get_subsimplex( n, 0, t, v1 );
            entry.column = mesh.get_subsimplex( n, 0, t, v2 );
            
            if( v1 == v2 )
                entry.value = 2. * factorial_numerical(n) * measure / factorial_numerical( 2 + n );
            else
                entry.value =      factorial_numerical(n) * measure / factorial_numerical( 2 + n );
                    
            if( mesh.get_flag( 0, vertex1 ) == SimplexFlag::SimplexFlagDirichlet or mesh.get_flag( 0, vertex2 ) == SimplexFlag::SimplexFlagDirichlet )
                entry.value = 0.;

            ret.setentry( index_of_entry, entry );
            
        }
            
    }

    ret.sortandcompressentries();

    return ret;
}




SparseMatrix LagrangeBrokenStiffnessMatrix( const Mesh& mesh, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_volumes = mesh.count_simplices( n );
    
    const int dim_in  = num_volumes * (n+1);
    const int dim_out = num_volumes * (n+1);
    
    
    // Set up sparse matrix
    
    SparseMatrix ret( dim_out, dim_in, num_volumes * (n+1)*(n+1) );
    
    
    // go over the n simplices and their k subsimplices
    for( int t  = 0; t  <  num_volumes; t++  )
    {
        
        DenseMatrix Jac = mesh.getTransformationJacobian( n, t );
        
        DenseMatrix GradProds = mesh.getGradientProductMatrix( n, t );
        
        Float measure = mesh.getMeasure( n, t );
        
        for( int v1 = 0; v1 <= n;           v1++ )
        for( int v2 = 0; v2 <= n;           v2++ )
        {
            
            SparseMatrix::MatrixEntry entry;
            
            int index_of_entry = t * (n+1)*(n+1) + v1 * (n+1) + v2;
            
            entry.row    = t * (n+1) + v1;
            
            entry.column = t * (n+1) + v2;
            
            // entry.value  = 0.;
            
            entry.value = GradProds( v1, v2 ) * measure; // /factorial_numerical( n );
            
            ret.setentry( index_of_entry, entry );
            
        }
        
    }
    
    return ret;
}





SparseMatrix LagrangeStiffnessMatrix( const Mesh& mesh, int r )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_volumes = mesh.count_simplices( n );
        
    const int num_vertices = mesh.count_simplices( 0 );
    
    const int dim_in  = num_vertices;
    const int dim_out = num_vertices;
    
    
    // Set up sparse matrix
    
    SparseMatrix ret( dim_out, dim_in, num_volumes * (n+1)*(n+1) );
    
    
    // go over the n simplices and their k subsimplices
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int t  = 0; t  <  num_volumes; t++  )
    for( int v1 = 0; v1 <= n;           v1++ )
    for( int v2 = 0; v2 <= n;           v2++ )
    {
        
        SparseMatrix::MatrixEntry entry;
        
        int index_of_entry = t * (n+1)*(n+1) + v1 * (n+1) + v2;
        
        int vertex1 = mesh.get_subsimplex( n, 0, t, v1 );
        int vertex2 = mesh.get_subsimplex( n, 0, t, v2 );

        entry.row    = vertex1; 
        entry.column = vertex2;
        
        //DenseMatrix Jac = mesh.getTransformationJacobian( n, t );
        
        DenseMatrix GradProds = mesh.getGradientProductMatrix( n, t );
        
        Float measure = mesh.getMeasure( n, t );
        
        entry.value = GradProds( v1, v2 ) * measure; // / factorial_numerical( n );
        
        if( mesh.get_flag( 0, vertex1 ) == SimplexFlag::SimplexFlagDirichlet or mesh.get_flag( 0, vertex2 ) == SimplexFlag::SimplexFlagDirichlet )
            entry.value = 0.;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}





SparseMatrix LagrangeInclusionMatrix( const Mesh& mesh, int n, int r )
{
    
    // check whether the parameters are right 
    assert( n >= 0 && n <= mesh.getinnerdimension() );
    assert( r >= 1 );
    assert( r == 1 );
    
    const int num_simplices = mesh.count_simplices( n );
    const int num_vertices  = mesh.count_simplices( 0 );
    
    const int dim_out = num_simplices * (n+1);
    const int dim_in  = num_vertices;
    
    const int num_entries = num_simplices * (n+1);
    
    

    SparseMatrix ret( dim_out, dim_in, num_entries );
    
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s  = 0; s  <  num_simplices; s++  )
    for( int vi = 0; vi <= n;             vi++ )
    {
        
        int index_of_entry = s * (n+1) + vi; 
        
        int vertex = mesh.get_subsimplex( n, 0, s, vi );
        
        SparseMatrix::MatrixEntry entry;
        
        entry.row    = s * (n+1) + vi;
        entry.column = vertex;
        entry.value  = 1.0;
        
        if( mesh.get_flag( 0, vertex ) == SimplexFlag::SimplexFlagDirichlet )
            entry.value = 0.0;
        
        ret.setentry( index_of_entry, entry );
        
    }
    
    return ret;
}





MatrixCSR LagrangeCoefficientMassMatrix( const Mesh& mesh, int r, int w, const std::function<Float(const FloatVector&)> weight )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_vertices = mesh.count_simplices( 0 );
    
    const int num_edges    = mesh.count_simplices( 1 );
    
    const int num_volumes  = mesh.count_simplices( n );
        
    const int dim_in  = num_vertices;
    const int dim_out = num_vertices;
    
    
    LOG << "Setting up A..." << num_vertices << space << nl;

    std::vector<int> A( num_vertices + 1 );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int e = 0; e < num_edges; e++ )
    {
        int v1 = mesh.get_subsimplex(1,0,e,0);
        int v2 = mesh.get_subsimplex(1,0,e,1);
        
        assert( 0 <= v1 && v1 < num_vertices );
        assert( 0 <= v2 && v2 < num_vertices );

        #if defined(_OPENMP)
        #pragma omp atomic
        #endif
        A[v1+1]++;

        #if defined(_OPENMP)
        #pragma omp atomic
        #endif
        A[v2+1]++;
    }   

    for( int v = 1; v <= num_vertices; v++ ) A[v] = A[v] + 1;

    LOG << "max width... ";
    int max_width = 0;
    for( int v = 1; v <= num_vertices; v++ ) max_width = maximum( max_width, A[v] );
    LOG << max_width << nl;

    A[0] = 0; 
    for( int v = 1; v <= num_vertices; v++ )
        A[v] += A[v-1];

    const int num_entries = A.back();

    for( int v = 1; v <= num_vertices; v++ ) assert( A[v] > A[v-1] );

    
    // set up the column indices     

    LOG << "Setting up C..." << num_entries << nl;

    std::vector<int> C( num_entries );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int v = 0; v < num_vertices; v++ )
    {
        int counter = 0;
        
        C[ A[v] + (counter++) ] = v;
        
        auto parent_edges = mesh.getsupersimplices(1,0,v);

        for( const auto e : parent_edges )
        {
            
            int v1 = mesh.get_subsimplex(1,0,e,0);
            if( v1 != v ) C[ A[v] + (counter++) ] = v1;
            
            int v2 = mesh.get_subsimplex(1,0,e,1);
            if( v2 != v ) C[ A[v] + (counter++) ] = v2;

            assert( v == v1 or v == v2 );

        }

        assert( counter == A[v+1] - A[v] );

        // std::sort( A.begin() + A[v], A.begin() + A[v+1] );
        for( int i = A[v]+1; i < A[v+1]; i++ )
        for( int j = A[v]+1; j < A[v+1]; j++ )
            if( C[j-1] > C[j] )
                std::swap( C[j-1], C[j] );

        for( int i = A[v]+1; i < A[v+1]; i++ )
            Assert( C[i] > C[i-1], C[i], C[i-1], counter, A[v+1], A[v] );
                
    }

    
    LOG << "Setting up V..." << num_entries << nl;

    std::vector<Float> V( num_entries, 0. );


    // assemble algebraic auxiliary material
    // - lagrange points in barycentric coordinates 
    // - coefficients of Lagrange polynomials
    // mass matrices 

    const auto lpbcs = InterpolationPointsInBarycentricCoordinates( n, w );

    const auto lpcoeff = Inverse( PointValuesOfMonomials( w, lpbcs ) );
    
    const auto polymassmatrix_per_point = polynomialmassmatrices_per_lagrangepoint( n, r, w );
    

    // go over the n simplices and their vertices
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_volumes; s++ )
    {

        // assemble some data for the element 
        // - measure 
        // - barycentric coordinates 
        // - lagrange points 

        Float measure = mesh.getMeasure( n, s );

        const auto vertex_coordinates = mesh.getVertexCoordinateMatrix( n, s );
        const auto lpeucl             = vertex_coordinates * lpbcs;

        DenseMatrix local_mass_matrix( n+1, n+1, 0. );

        for( int p = 0; p < polymassmatrix_per_point.size(); p++ )
        {
            const auto& polyMM = polymassmatrix_per_point[p];
            
            const Float weight_at_point = weight( lpeucl.getcolumn(p) );

            if( w == 0 ) assert( ( polyMM - polynomialmassmatrix(n,r) ).is_numerically_small() );

            auto fullMM = polyMM;
            fullMM *= ( weight_at_point * measure );

            local_mass_matrix += fullMM;
        }

        assert( local_mass_matrix.isfinite() );

        for( int v1 = 0; v1 <= n; v1++ )
        for( int v2 = 0; v2 <= n; v2++ )
        {

            int vertex1 = mesh.get_subsimplex( n, 0, s, v1 );
            int vertex2 = mesh.get_subsimplex( n, 0, s, v2 );

            if( mesh.get_flag( 0, vertex1 ) == SimplexFlag::SimplexFlagDirichlet or mesh.get_flag( 0, vertex2 ) == SimplexFlag::SimplexFlagDirichlet )
                local_mass_matrix( v1, v2 ) = 0.;

        }

        assert( local_mass_matrix.isfinite() );

        for( int v1 = 0; v1 <= n; v1++ )
        for( int v2 = 0; v2 <= n; v2++ )
        {
            int vertex1 = mesh.get_subsimplex( n, 0, s, v1 );
            int vertex2 = mesh.get_subsimplex( n, 0, s, v2 );

            Assert( 0 <= vertex1 && vertex1 < num_entries, vertex1 );
            Assert( 0 <= vertex2 && vertex2 < num_entries, vertex2 );

            for( int i = A[vertex1]; i < A[vertex1+1]; i++ ) {
            
                Assert( 0 <= i && i < num_entries, i );
            
                if( C[i] == vertex2 ) {
                    
                    #if defined(_OPENMP)
                    #pragma omp atomic
                    #endif
                    V[i] += local_mass_matrix( v1, v2 );

                    assert( std::isfinite( V[i] ) );

                }

            }
                    
        }

    }

    // the data are complete

    for( const auto& v : V ) assert( std::isfinite(v) );

    return MatrixCSR( dim_out, dim_in, std::move(A), std::move(C), std::move(V) );

}







MatrixCSR LagrangeCoefficientStiffnessMatrix( const Mesh& mesh, int r, int w, const std::function<DenseMatrix(const FloatVector&)> weight )
{
    
    // check whether the parameters are right 
    // only lowest order here
    
    assert( r >= 0 );
    assert( r == 1 ); // only lowest order for the time being
    
    // Auxiliary calculations and preparations
    
    int n = mesh.getinnerdimension();
    
    const int num_volumes  = mesh.count_simplices( n );
        
    const int num_edges    = mesh.count_simplices( 1 );
        
    const int num_vertices = mesh.count_simplices( 0 );
    
    const int dim_in  = num_vertices;
    const int dim_out = num_vertices;

    // Compute the number of entries
    // and set up A

    LOG << "Setting up A..." << num_vertices << nl;

    std::vector<int> A( num_vertices + 1 );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int e = 0; e < num_edges; e++ )
    {
        int v1 = mesh.get_subsimplex(1,0,e,0);
        int v2 = mesh.get_subsimplex(1,0,e,1);
        
        assert( 0 <= v1 && v1 < num_vertices );
        assert( 0 <= v2 && v2 < num_vertices );

        #if defined(_OPENMP)
        #pragma omp atomic
        #endif
        A[v1+1]++;

        #if defined(_OPENMP)
        #pragma omp atomic
        #endif
        A[v2+1]++;
    }   

    for( int v = 1; v <= num_vertices; v++ ) A[v] = A[v] + 1;

    LOG << "max width... ";
    int max_width = 0;
    for( int v = 1; v <= num_vertices; v++ ) max_width = maximum( max_width, A[v] );
    LOG << max_width << nl;

    
    A[0] = 0; 
    for( int v = 1; v <= num_vertices; v++ )
        A[v] += A[v-1];

    const int num_entries = A.back();

    for( int v = 1; v <= num_vertices; v++ ) assert( A[v] > A[v-1] );

    
    // set up the column indices     

    LOG << "Setting up C..." << num_entries << nl;

    std::vector<int> C( num_entries );

    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int v = 0; v < num_vertices; v++ )
    {
        int counter = 0;
        
        C[ A[v] + (counter++) ] = v;
        
        auto parent_edges = mesh.getsupersimplices(1,0,v);

        for( const auto e : parent_edges )
        {
            
            int v1 = mesh.get_subsimplex(1,0,e,0);
            if( v1 != v ) C[ A[v] + (counter++) ] = v1;
            
            int v2 = mesh.get_subsimplex(1,0,e,1);
            if( v2 != v ) C[ A[v] + (counter++) ] = v2;

            assert( v == v1 or v == v2 );

        }

        assert( counter == A[v+1] - A[v] );

        // std::sort( A.begin() + A[v], A.begin() + A[v+1] );
        for( int i = A[v]+1; i < A[v+1]; i++ )
        for( int j = A[v]+1; j < A[v+1]; j++ )
            if( C[j-1] > C[j] )
                std::swap( C[j-1], C[j] );

        for( int i = A[v]+1; i < A[v+1]; i++ )
            Assert( C[i] > C[i-1], C[i], C[i-1], counter, A[v+1], A[v] );
                
    }

    
    
    LOG << "Setting up V..." << num_entries << nl;

    std::vector<Float> V( num_entries, 0. );

    // assemble algebraic auxiliary material
    // - lagrange points in barycentric coordinates 
    // - coefficients of Lagrange polynomials
    // mass matrices 

    const auto lpbcs = InterpolationPointsInBarycentricCoordinates( n, w );

    const auto lpcoeff = Inverse( PointValuesOfMonomials( w, lpbcs ) );
    
    const auto polymassmatrix_per_point = polynomialmassmatrices_per_lagrangepoint( n, r, w );
    

    // go over the n simplices and their vertices
    #if defined(_OPENMP)
    #pragma omp parallel for
    #endif
    for( int s = 0; s < num_volumes; s++ )
    {

        // assemble some data for the element 
        // - measure 
        // - barycentric coordinates 
        // - lagrange points 

        Float measure = mesh.getMeasure( n, s );
        assert( measure >= 0. );

        const DenseMatrix  GM    = mesh.getGradientMatrix( n, s );
        const DenseMatrix& extGM = GM; //SubdeterminantMatrix( GM, 1 );

        const auto vertex_coordinates = mesh.getVertexCoordinateMatrix( n, s );
        const auto lpeucl             = vertex_coordinates * lpbcs;

        DenseMatrix local_mass_matrix( n+1, n+1, 0. );





        // compute the mass matrix contribution 
        // for each lagrange point 
        
        for( int p = 0; p < polymassmatrix_per_point.size(); p++ )
        {
            const auto& polyMM = polymassmatrix_per_point[p];
            
            const DenseMatrix matrix_at_point = weight( lpeucl.getcolumn(p) );

            const auto formMM = Transpose(extGM) * matrix_at_point * extGM;
            // const auto formMM = MatrixTripleMult( matrix_at_point, extGM );

            // DenseMatrix GPM = SubdeterminantMatrix( mesh.getGradientProductMatrix( n, s ), k );
            // assert( ( GPM - formMM ).is_numerically_small() );

            if( w == 0 ) assert( ( polyMM - polynomialmassmatrix(n,r) ).is_numerically_small() );

            //auto fullMM = measure * MatrixTensorProduct( polyMM, formMM );
            auto fullMM = MatrixTensorProduct( polyMM, formMM );
            fullMM *= measure;

            local_mass_matrix += fullMM;
        }




        assert( local_mass_matrix.isfinite() );

        for( int v1 = 0; v1 <= n; v1++ )
        for( int v2 = 0; v2 <= n; v2++ )
        {

            int vertex1 = mesh.get_subsimplex( n, 0, s, v1 );
            int vertex2 = mesh.get_subsimplex( n, 0, s, v2 );

            if( mesh.get_flag( 0, vertex1 ) == SimplexFlag::SimplexFlagDirichlet or mesh.get_flag( 0, vertex2 ) == SimplexFlag::SimplexFlagDirichlet )
                local_mass_matrix( v1, v2 ) = 0.;

        }

        assert( local_mass_matrix.isfinite() );

        for( int v1 = 0; v1 <= n; v1++ )
        for( int v2 = 0; v2 <= n; v2++ )
        {
            int vertex1 = mesh.get_subsimplex( n, 0, s, v1 );
            int vertex2 = mesh.get_subsimplex( n, 0, s, v2 );

            Assert( 0 <= vertex1 && vertex1 < num_entries, vertex1 );
            Assert( 0 <= vertex2 && vertex2 < num_entries, vertex2 );

            for( int i = A[vertex1]; i < A[vertex1+1]; i++ ) {

                Assert( 0 <= i && i < num_entries, i );
                
                if( C[i] == vertex2 ) {
                    
                    #if defined(_OPENMP)
                    #pragma omp atomic
                    #endif
                    V[i] += local_mass_matrix( v1, v2 );

                    assert( std::isfinite( V[i] ) );
                    
                }
            
            }
                    
        }

    }

    // the data are complete

    for( const auto& v : V ) assert( std::isfinite(v) );

    return MatrixCSR( dim_out, dim_in, std::move(A), std::move(C), std::move(V) );

}

