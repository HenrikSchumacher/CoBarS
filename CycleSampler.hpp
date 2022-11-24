#ifndef CYCLE_SAMPLER_HPP
    #define CYCLE_SAMPLER_HPP

    #include <cmath>
    #include <algorithm>
    #include <iostream>
    #include <random>
    #include <cstring>
    #include <array>


    #define EIGEN_NO_DEBUG
    #define EIGEN_USE_BLAS
    #define EIGEN_USE_LAPACKE
    #include <eigen3/Eigen/Dense>

    #include "Tensors/Tensors.hpp"

<<<<<<< HEAD
    namespace CycleSampler {
=======
    namespace CycleSampler
    {
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
        
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

    #include "src/ShiftMap.hpp"
    #include "src/SamplerBase.hpp"
    #include "src/Sampler.hpp"

    #include "src/RandomVariables.hpp"


    #include "src/MomentPolytopeSampler.hpp"

#endif
