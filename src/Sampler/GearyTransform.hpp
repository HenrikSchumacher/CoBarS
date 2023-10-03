private:

    struct GearyTransform
    {
        const Real mean_X;
        const Real mean_Y;
        const Real var_X;
        const Real cov_XY;
        const Real var_Y;
        
        GearyTransform( const Real mean_X_, const Real mean_Y_,
                        const Real var_X_,  const Real cov_XY_, const Real var_Y_ )
        :   mean_X ( mean_X_ )
        ,   mean_Y ( mean_Y_ )
        ,   var_X  ( var_X_  )
        ,   cov_XY ( cov_XY_ )
        ,   var_Y  ( var_Y_  )
        {}
        
        Real operator()( const Real T ) const
        {
            const Real numerator = mean_Y * T - mean_X;
            const Real denominator_squared = Abs( var_Y * T * T + var_X - two * cov_XY * T );
            const Real denominator = Sqrt( denominator_squared );
            
            return numerator / denominator;
        }
        
        Real D( const Real T ) const
        {
            const Real numerator = mean_Y * T - mean_X;
            const Real denominator_squared = Abs( var_Y * T * T + var_X - two * cov_XY * T );
            const Real denominator = Sqrt( denominator_squared );
            
            return ( mean_Y * denominator - numerator * two * (var_Y * T - cov_XY) ) / denominator_squared;
        }
    };
