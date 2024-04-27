#ifndef INCLUDEGUARD_FEM_FINITEDIFF_HPP
#define INCLUDEGUARD_FEM_FINITEDIFF_HPP


#include <functional>
#include <vector>

#include "../basic.hpp"
#include "../utility/stl.hpp"
#include "../combinatorics/indexrange.hpp"
#include "../combinatorics/indexmap.hpp"
#include "../combinatorics/multiindex.hpp"
#include "../combinatorics/generatemultiindices.hpp"
#include "../combinatorics/generateindexmaps.hpp"
#include "../operators/floatvector.hpp"
#include "../dense/densematrix.hpp"



class AlternatingForm
{
    
    private:
        int d;
        int k;
        std::function<FloatVector(const FloatVector&)> formula;   
    
    public:
    
        AlternatingForm( int d, int k, const std::function<FloatVector(const FloatVector&)>& formula )
        : d(d), k(k), formula(formula)
        {
            check();
        }   
        
        AlternatingForm( const AlternatingForm& ) = default;
        AlternatingForm( AlternatingForm&& ) = default;

        virtual ~AlternatingForm() = default;
        
        void check()
        {
            assert( d >= 0 && 0 <= k && k <= d );
        }
        
        
        int getdimension() const
        {
            return d;
        }
        
        int getdegree() const
        {
            return k;
        }
        
        FloatVector operator()( const FloatVector& point )
        {
            assert( point.getdimension() == getdimension() );
            FloatVector ret = formula( point );
            assert( ret.getdimension() == binomial_integer(d,k) );
            return ret;
        }
        
        
        AlternatingForm exteriorderivative( Float stepsize = desired_closeness ) const 
        {
            if( stepsize < 0 ) stepsize = std::cbrt( std::numeric_limits<Float>::epsilon() );

            const std::vector<IndexMap> sigmas_dst = generateSigmas( IndexRange( 1, k+1 ), IndexRange( 0, d-1 ) );
            const std::vector<IndexMap> sigmas_src  = generateSigmas( IndexRange( 1, k   ), IndexRange( 0, d-1 ) );
            
            std::vector<DenseMatrix> pattern( d, DenseMatrix( sigmas_dst.size(), sigmas_src.size(), 0 ) );
            
            assert( sigmas_dst.size() == binomial_integer( d, k+1 ) );
            assert( sigmas_src.size()  == binomial_integer( d, k   ) );
            
            for( int src_form_index = 0; src_form_index < sigmas_src.size(); src_form_index++ )
            for( int p = 0; p < d; p++ ) {
                
                const IndexMap& src_form = sigmas_src[src_form_index];
                
                assert( pattern[p].getdimin()  == sigmas_src.size() );
                assert( pattern[p].getdimout() == sigmas_dst.size() );
                
                if( not src_form.has_value_in_range(p) ) {
                    
                    IndexMap new_form = expand_one( src_form, p );
                    int new_form_index = find_index( sigmas_dst, new_form );
                    pattern.at(p).at( new_form_index, src_form_index ) = sign_power( new_form.preimageof(p) );
                    
                }
                
            }
            
            int d = this->d;
            int k = this->k;
            const auto formula = this->formula;
            
            // const auto newformula = [k,d,pattern,formula,stepsize]( const FloatVector& point ) -> FloatVector {
            const auto newformula = [=]( const FloatVector& point ) -> FloatVector {
                
                int dim_src = binomial_integer(d,k  );
                int dim_dst = binomial_integer(d,k+1);
                
                FloatVector ret( dim_dst, 0. );
                
//                 LOG << dim_src << space << dim_dst << space << pattern.size() << nl;
                
                assert( point.getdimension() == d );
                assert( pattern.size() == d );
                
                std::vector<FloatVector> diffs; diffs.reserve(d);
                
                for( int i = 0; i < d; i++ )
                {
                    FloatVector forward  = formula( point + stepsize * unitvector(d,i) );
                    FloatVector backward = formula( point - stepsize * unitvector(d,i) );
                    FloatVector diff = 0.5 * ( forward - backward ) / stepsize;
                    diffs.push_back( diff );
                    assert( diff.getdimension() == dim_src );
                }
                
                for( int i = 0; i < d; i++ )
                    ret += pattern[i] * diffs[i];
                    
                return ret;
                
            };
    
            return AlternatingForm( d, k+1, newformula );
            
        }
        
        
        
        AlternatingForm laplacian( Float stepsize = desired_closeness ) const 
        {
            
            if( stepsize < 0 ) stepsize = std::cbrt( std::numeric_limits<Float>::epsilon() );

            int d = this->d;
            int k = this->k;
            const auto formula = this->formula;
            
            const auto newformula = [d,k,formula,stepsize]( const FloatVector& point ) -> FloatVector {
                
                int fielddim = binomial_integer(d,k);
                
                FloatVector ret( fielddim, 0. );
                
                FloatVector mid = formula( point );
                
                std::vector<FloatVector> diffs; diffs.reserve(d);
                
                for( int i = 0; i < d; i++ )
                {
                    FloatVector forward  = formula( point + stepsize * unitvector(d,i) );
                    FloatVector backward = formula( point - stepsize * unitvector(d,i) );
                    diffs.push_back( 2 * mid - forward - backward );
                }
                
                for( int i = 0; i < d; i++ ) ret = ret - diffs[i];
                
                return ret / ( stepsize * stepsize );
                
            };
    
            return AlternatingForm( d, k, newformula );
            
        }
        
        
        
};








#endif
