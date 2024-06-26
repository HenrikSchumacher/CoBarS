#pragma once

namespace CoBarS
{
    template<typename Sampler_T> class RandomVariable;

    
    /*!
     * @brief A struct to carry the options for `CoBarS::SamplerBase` and `CoBarS::Sampler`.
     */
    
    template<typename Real, typename Int>
    struct SamplerSettings
    {
//        Real tolerance            = std::sqrt(Scalar::eps<Real>);
        Real tolerance            = std::pow(Scalar::eps<Real>,Frac<Real>(3,4));
        Real give_up_tolerance    = 128 * Scalar::eps<Real>;
        Real regularization       = static_cast<Real>(0.01);
        Int  max_iter             = 1000;
        
        Real Armijo_slope_factor  = static_cast<Real>(0.01);
        Real Armijo_shrink_factor = static_cast<Real>(0.5);
        Int  max_backtrackings    = 20;
        
        bool use_linesearch       = true;
        
        SamplerSettings() {}
        
        ~SamplerSettings() = default;
        
        SamplerSettings( const SamplerSettings & other )
        :   tolerance(other.tolerance)
        ,   give_up_tolerance(other.give_up_tolerance)
        ,   regularization(other.regularization)
        ,   max_iter(other.max_iter)
        ,   Armijo_slope_factor(other.Armijo_slope_factor)
        ,   Armijo_shrink_factor(other.Armijo_shrink_factor)
        ,   max_backtrackings(other.max_backtrackings)
        ,   use_linesearch(other.use_linesearch)
        {}
        
        void PrintStats() const
        {
            valprint( "tolerance           ", tolerance           , 16 );
            valprint( "give_up_tolerance   ", give_up_tolerance   , 16 );
            valprint( "regularization      ", regularization      , 16 );
            valprint( "max_iter            ", max_iter            , 16 );
            valprint( "Armijo_slope_factor ", Armijo_slope_factor , 16 );
            valprint( "Armijo_shrink_factor", Armijo_shrink_factor, 16 );
            valprint( "max_backtrackings   ", max_backtrackings   , 16 );
            valprint( "use_linesearch      ", use_linesearch      , 16 );
        }
    };
    
    
} // namespace CoBarS
