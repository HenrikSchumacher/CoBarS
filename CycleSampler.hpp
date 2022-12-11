#ifndef CYCLE_SAMPLER_HPP
    #define CYCLE_SAMPLER_HPP

    #include <cmath>
    #include <algorithm>
    #include <iostream>
    #include <random>
    #include <cstring>
    #include <array>

    #include "Tensors/Tensors.hpp"

    namespace CycleSampler {
        
        using namespace Tools;
        using namespace Tensors;
        
    } // namespace CycleSampler

    //extern "C"
    //{
    //    #include <gsl/gsl_rng.h>
    //    #include <gsl/gsl_randist.h>
    //}
    //
    //#include <plcTopology.h>
    //
    //int PD_VERBOSE = 0;


    #include "src/Sampler.hpp"

    #include "src/RandomVariable.hpp"

    #include "src/RandomVariables/BarycenterNorm.hpp"   // Should always evaluate to 0.
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
