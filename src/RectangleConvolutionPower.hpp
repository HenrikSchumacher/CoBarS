#pragma once

namespace CycleSampler
{

    template<typename Real_,typename Int_>
    class RectangleConvolutionPower
    {
        ASSERT_REAL(Real_);
        ASSERT_INT(Int_);
        
    public:
        
        using Real = Real_;
        using Int = Int_;
        
        RectangleConvolutionPower() = default;
        
        explicit RectangleConvolutionPower( Int n_ )
        :   n         ( n_                   )
        ,   n_Real    ( static_cast<Real>(n) )
        ,   a         ( n_                   )
        {}
        
        ~RectangleConvolutionPower() = default;
        
        
        static constexpr Real factor = cSqrt( Scalar::Pi<Real> * Scalar::Half<Real>);
        
        Real Value( const Real x )
        {
            if( (x <= - n_Real) || (x >= n_Real) )
            {
                return Scalar::Zero<Real>;
            }
            
            const Real x_plus_n = x + n_Real;
            
            const Int i_0 = static_cast<Int>(std::floor(Scalar::Half<Real> * (x_plus_n)));
            
            a.SetZero();
            
            a[i_0] = Scalar::One<Real>;
            
            for( Int k = 1; k < n; ++k )
            {
                const Int i_begin = std::max( Scalar::Zero<Int>, i_0 - k );
                const Int i_end   = std::min( i_0 + 1, n - k );
                
                const Real two_k = static_cast<Real>( 2 * k );
                
                const Real two_k_plus_2 = two_k + Scalar::Two<Real>;
                
                const Real scale = Scalar::Inv<Real>( two_k );
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    const Real t = x_plus_n - static_cast<Real>( 2 * i );
                    
                    a[i] = scale * ( a[i] * t + a[i+1] * ( two_k_plus_2 - t ) );
                }
            }
            
            return factor * a[0];
        }
        
        Real Derivative( const Real x )
        {
            if( (x <= - n_Real) || (x >= n_Real) )
            {
                return Scalar::Zero<Real>;
            }
            
            const Real x_plus_n = x + n_Real;
            
            const Int i_0 = static_cast<Int>(std::floor(Scalar::Half<Real> * (x_plus_n)));
            
            a.SetZero();
            
            a[i_0] = Scalar::Half<Real>; // <-- We don't have to multiply with 0.5 in the end.
            
            for( Int k = 1; k < n-1; ++k )
            {
                const Int i_begin = std::max( Scalar::Zero<Int>, i_0 - k );
                const Int i_end   = std::min( i_0 + 1, n - k );
                
                const Real two_k = static_cast<Real>( 2 * k );
                
                const Real two_k_plus_2 = two_k + Scalar::Two<Real>;
                
                const Real scale = Scalar::Inv<Real>( two_k );
                
                for( Int i = i_begin; i < i_end; ++i )
                {
                    const Real t = x_plus_n - static_cast<Real>( 2 * i );
                    
                    a[i] = scale * ( a[i] * t + a[i+1] * ( two_k_plus_2 - t ) );
                }
            }
            
            return factor * (a[0] - a[1]);
        }
        
        Int GetPower() const
        {
            return n;
        }
        
        Int SetPower( const Int n_ )
        {
            if( n_ > 0 )
            {
                n = n_;
                n_Real = scalar_cast<Real>(n_);
                a.Resize(n);
            }
        }
        
    private:
        
        Int  n       = 0;
        Real n_Real  = 0;
        
        Tensor1<Real,Int> a;
        
    }; // class RectangleConvolutionPower
    
} // namespace CycleSampler
