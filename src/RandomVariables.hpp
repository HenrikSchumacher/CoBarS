#pragma once


#include "RandomVariables/RandomVariableBase.hpp"
#include "RandomVariables/RandomVariable.hpp"

namespace CyclicSampler
{
    #include "RandomVariables/BarycenterNorm.hpp"   // Should always evaluate to 0.
    #include "RandomVariables/ChordLength.hpp"
    #include "RandomVariables/DiagonalLength.hpp"
    #include "RandomVariables/Gyradius.hpp"
    #include "RandomVariables/GyradiusP.hpp"
    #include "RandomVariables/HydrodynamicRadius.hpp"
    #include "RandomVariables/ShiftNorm.hpp"
    #include "RandomVariables/TotalCurvature.hpp"
    #include "RandomVariables/BendingEnergy.hpp"
    #include "RandomVariables/MaxAngle.hpp"
    #include "RandomVariables/EdgeSpaceSamplingWeight.hpp"
    #include "RandomVariables/EdgeQuotientSpaceSamplingWeight.hpp"
    #include "RandomVariables/IterationCount.hpp"
        
        
    // Add your own random variables here.
    #include "RandomVariables/ExampleFunction.hpp"
}
