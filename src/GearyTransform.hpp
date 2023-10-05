#pragma once

namespace CycleSampler
{
    
    template<typename Real>
    struct GearyTransform
    {
        const Real mean_X;
        const Real mean_Y;
        const Real var_X;
        const Real cov_X_Y;
        const Real var_Y;
        
        GearyTransform( const Real mean_X_, const Real mean_Y_,
                        const Real var_X_,  const Real cov_X_Y_, const Real var_Y_
        )
        :   mean_X ( mean_X_   )
        ,   mean_Y ( mean_Y_   )
        ,   var_X  ( var_X_    )
        ,   cov_X_Y ( cov_X_Y_ )
        ,   var_Y  ( var_Y_    )
        {}
        
        GearyTransform( GearyTransform & other )
        :   mean_X  ( other.mean_X  )
        ,   mean_Y  ( other.mean_Y  )
        ,   var_X   ( other.var_X   )
        ,   cov_X_Y ( other.cov_X_Y )
        ,   var_Y   ( other.var_Y   )
        {}
        
        Real operator()( const Real T ) const
        {
            const Real numerator = mean_Y * T - mean_X;
            const Real denominator_squared = std::abs( var_Y * T * T + var_X - Scalar::Two<Real> * cov_X_Y * T );
            const Real denominator = std::sqrt( denominator_squared );
            
            return numerator / denominator;
        }
        
        Real D( const Real T ) const
        {
            const Real numerator = mean_Y * T - mean_X;
            const Real denominator_squared = Abs( var_Y * T * T + var_X - Scalar::Two<Real> * cov_X_Y * T );
            const Real denominator = Sqrt( denominator_squared );
            
            return ( mean_Y * denominator - numerator * Scalar::Two<Real> * (var_Y * T - cov_X_Y) ) / denominator_squared;
        }
    };
    
} // namespace CycleSampler
