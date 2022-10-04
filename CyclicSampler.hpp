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

//extern "C"
//{
//    #include <gsl/gsl_rng.h>
//    #include <gsl/gsl_randist.h>
//}
//
//#include <plcTopology.h>
//
//int PD_VERBOSE = 0;

#include "src/SmallSymmetricMatrix.hpp"

#include "src/ShiftMap.hpp"
#include "src/CyclicSamplerBase.hpp"
#include "src/CyclicSampler.hpp"

#include "src/RandomVariables.hpp"


#include "src/MomentPolytopeSampler.hpp"
