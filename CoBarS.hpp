#ifndef COBARS_HPP

    #define COBARS_HPP


/*! \mainpage My Personal Index Page
 *
 * \section intro_sec Introduction
 *
 * This is the introduction.
 *
 * \section install_sec Installation
 *
 * \subsection step1 Step 1: Opening the box
 *
 * etc...
 */

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
        
/*!
 * Commulative distribution function (CDF) of normal distribution.
 *
 * @tparam Real A real floating point type.
 *
 * @param z The argument.
 */
        
        template<typename Real>
        inline Real N_CDF( const Real z )
        {
            static_assert(FloatQ<Real>, "");
            
            constexpr Real threshold = static_cast<Real>(8);
            
            constexpr Real factor = Inv( cSqrt( Scalar::Two<Real> ) );
            
            return ( z < -threshold ) ? Scalar::Zero<Real> :  ( z > threshold ) ? Scalar::One<Real> : (Scalar::Half<Real> + Scalar::Half<Real> * std::erf( factor * z ));
        }
        
/*!
 * Probability distribution function (PDF) of normal distribution.
 *
 * @tparam Real A real floating point type.
 *
 * @param z The argument.
 */
        template<typename Real>
        inline Real N_PDF( const Real z )
        {
            static_assert(FloatQ<Real>, "");
            
            constexpr Real factor = Inv( cSqrt( Scalar::TwoPi<Real> ) );

            return factor * std::exp( - Scalar::Half<Real> * z * z );
        }
        
    } // namespace CoBarS
    
    #include "src/MT64.hpp"
    #include "src/Xoshiro256Plus.hpp"
    #include "src/PCG64.hpp"
    #include "src/WyRand.hpp"

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

    #include "src/RandomVariables/IterationCount.hpp"
        
        
    // Add your own random variables here!
    #include "src/RandomVariables/ExampleFunction.hpp"

    #include "src/ActionAngleSampler.hpp"

    #include "src/DouadyEarleExtension.hpp"


#endif
