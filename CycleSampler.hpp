#ifndef CYCLE_SAMPLER_HPP
    #define CYCLE_SAMPLER_HPP

    #include <cmath>
    #include <algorithm>
    #include <iostream>
    #include <random>
    #include <cstring>

    #include "Tensors/Tensors.hpp"

    namespace CycleSampler
    {
        
        using namespace Tools;
        using namespace Tensors;

#include "src/SplitMix64.hpp"
#include "src/MersenneTwister.hpp"
#include "src/Xoshiro256Plus.hpp"
    }



    #include "src/Sampler.hpp"

    #include "src/RandomVariable.hpp"

    #include "src/RandomVariables/BarycenterNorm.hpp"
    #include "src/RandomVariables/ChordLength.hpp"
    #include "src/RandomVariables/DiagonalLength.hpp"
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



    #include "src/MomentPolytopeSampler.hpp"

#endif
