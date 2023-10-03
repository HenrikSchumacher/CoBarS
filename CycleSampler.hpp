#ifndef CYCLE_SAMPLER_HPP
    #define CYCLE_SAMPLER_HPP

    #include <cmath>
    #include <algorithm>
    #include <iostream>
    #include <random>
    #include <cstring>

    #include "submodules/Tensors/Tensors.hpp"
    #include "submodules/erfinv/erfinv.c"

    #include <istream>
    #include <ostream>

    namespace CycleSampler
    {
        using namespace Tools;
        using namespace Tensors;
    }
    
    #include "src/MT64.hpp"
    #include "src/Xoshiro256Plus.hpp"
    #include "src/PCG64.hpp"

    #include "src/SamplerSettings.hpp"
    #include "src/SamplerBase.hpp"
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



    #include "src/ActionAngleSampler.hpp"

#endif
