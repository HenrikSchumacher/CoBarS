#ifndef CYCLE_SAMPLER_HPP
    #define CYCLE_SAMPLER_HPP

    #include <cmath>
    #include <algorithm>
    #include <iostream>
    #include <random>
    #include <cstring>

    #include "submodules/Tensors/Tensors.hpp"

    #include <istream>
    #include <ostream>

    namespace CoBarS
    {
        using namespace Tools;
        using namespace Tensors;
        
        
        template<typename Real>
        inline Real N_CDF( const Real z_ )
        {
            // CDF of normal distribution

            
            constexpr Real threshold = static_cast<Real>(8);
            
            constexpr Real factor = Inv( cSqrt( Scalar::Two<Real> ) );
            
            return ( z_ < -threshold ) ? Scalar::Zero<Real> :  ( z_ > threshold ) ? Scalar::One<Real> : (Scalar::Half<Real> + Scalar::Half<Real> * std::erf( factor * z_ ));
            
    //        return Scalar::Half<Real> * std::erf( factor * z_ );
        }
        
        template<typename Real>
        inline Real N_PDF( const Real z_ )
        {
            // PDF of normal distribution
            
            constexpr Real factor = Inv( cSqrt( Scalar::TwoPi<Real> ) );

            return factor * std::exp( - Scalar::Half<Real> * z_ * z_ );
        }
        
    } // namespace CoBarS
    
    #include "src/MT64.hpp"
    #include "src/Xoshiro256Plus.hpp"
    #include "src/PCG64.hpp"

    #include "src/GearyTransform.hpp"

    #include "src/SamplerSettings.hpp"
    #include "src/SamplerBase.hpp"
    #include "src/Sampler.hpp"

    #include "src/RandomVariable.hpp"

    #include "src/RandomVariables/BarycenterNorm.hpp"
    #include "src/RandomVariables/ChordLength.hpp"
    #include "src/RandomVariables/DiagonalLength.hpp"
    #include "src/RandomVariables/SquaredGyradius.hpp"
    #include "src/RandomVariables/Gyradius.hpp"
    #include "src/RandomVariables/GyradiusP.hpp"
    #include "src/RandomVariables/HydrodynamicRadius.hpp"
    #include "src/RandomVariables/ShiftNorm.hpp"
    #include "src/RandomVariables/TotalCurvature.hpp"
    #include "src/RandomVariables/BendingEnergy.hpp"
    #include "src/RandomVariables/MaxAngle.hpp"
    #include "src/RandomVariables/EdgeSpaceSamplingWeight.hpp"
    #include "src/RandomVariables/EdgeQuotientSpaceSamplingWeight.hpp"
    #include "src/RandomVariables/EdgeQuotientSpaceSamplingCorrection.hpp"

    #include "src/RandomVariables/IterationCount.hpp"
        
        
    // Add your own random variables here!
    #include "src/RandomVariables/ExampleFunction.hpp"



    #include "src/ActionAngleSampler.hpp"




#endif
