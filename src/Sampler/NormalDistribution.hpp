private:



    Real N_CDF( const Real z_ ) const
    {
        // CDF of normal distribution
        
        constexpr Real threshold = static_cast<Real>(8);
        
        constexpr Real factor = Inv( cSqrt( Scalar::Two<Real> ) );

        return ( z_ < -threshold ) ? Scalar::Zero<Real> :  ( z_ > threshold ) ? Scalar::One<Real> : Scalar::Half<Real> * std::erf( factor * z_ );
        
//        return Scalar::Half<Real> * std::erf( factor * z_ );
    }

    Real N_PDF( const Real z_ ) const
    {
        // PDF of normal distribution
        
        constexpr Real factor = Inv( cSqrt( Scalar::TwoPi<Real> ) );

        return factor * std::exp( - Scalar::Half<Real> * z_ * z_ );
    }
