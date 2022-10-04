#pragma once

#include <cmath>
#include <algorithm>
#include <iostream>
#include <random>
#include <cstring>
#include <array>



#define EIGEN_USE_BLAS
#define EIGEN_USE_LAPACKE
#include <eigen3/Eigen/Dense>

#include "Tensors/Tensors.hpp"

namespace CyclicSampler {
    using namespace Tools;
    using namespace Tensors;
} // namespace CyclicSampler

extern "C"
{
    #include <gsl/gsl_rng.h>
    #include <gsl/gsl_randist.h>
}

#include <plcTopology.h>

int PD_VERBOSE = 0;

#include "src/ShiftMap.hpp"
#include "src/CyclicSamplerBase.hpp"
#include "src/CyclicSampler.hpp"

#include "src/RandomVariables/RandomVariableBase.hpp"
#include "src/RandomVariables/RandomVariable.hpp"
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

// Only for debugging purposes:
#include "src/RandomVariables/BarycenterNorm.hpp"

#include "src/MomentPolytopeSampler.hpp"
