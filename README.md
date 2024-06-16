# CoBarS - Conformal Barycenter Sampling

by Jason Cantarella and Henrik Schumacher

A header-only C++ library for Monte-Carlo sampling of cylic polygons with prescribed edge lengths.

# Installation

TO make sure that all submodules are cloned, too, please clone by running the following in the command line:

    git clone --depth 1 --recurse-submodules --shallow-submodules git@github.com:HenrikSchumacher/CoBarS.git
        
    
To build the documentation, please install doxygen and then run

    doxygen Doxyfile
    
The documentation is then available in `doc/index.html`.
    
# Usage

_CoBarS_ is a header-only library with no dependencies other than the _C++ Standard Library_ and the _Standard Template Library_, and those that are already inlined to the repository (namely the pseudorandom number generators [PCG](https://github.com/imneme/pcg-cpp), [wy](https://github.com/alainesp/wy), and [Xoshiro256+](https://github.com/Reputeless/Xoshiro-cpp)).

Just include the header `CoBarS.hpp` into your C++ program via

    #include "CoBarS.hpp"    
        

You can find a more detailed version of the following example in the example program in `Example_RandomClosedPolygon/main.cpp`.
Suppose we want to create random, equilateral 64-gons along with their sampling weights. 
We first create a `CoBarS::Sampler` object; this is the main class of the package.
    
    constexpr int d = 3; // Dimensions of the ambient space has to be a compile-time constant.
    const int edge_count = 64;
    
    // Create an instance of the CoBarS sampler.
    // This assumes that equilateral polygons shall be created.
    CoBarS::Sampler<d,double,std::size_t> S (edge_count);

    
Next we need some buffers to store the random polygons and the sampling weights. You can allocated them yourself, but we use the container `std::vector` from the _STL_.
    
    const std::size_t sample_count      = 10000000;

    // Create containers for the polygons `p` and sampling weights `K`.
    std::vector<double> p ( sample_count * (edge_count + 1) * d );
    std::vector<double> K ( sample_count );

    
Finally submit the buffers to the member function `CreateRandomClosedPolygons`:

    // Whether we want the sampling weights for the quotient space of polygons modulo SO(d) (`quot_space_Q = true`) or not (`quot_space_Q = false`).
    const bool quot_space_Q = true;
    
    const std::size_t thread_count = 8;
    
    S.CreateRandomClosedPolygons(
        &p[0], &K[0], sample_count, quot_space_Q, thread_count
    );
    


# Compilation


We developped and tested _CoBarS_ most thoroughly with the _Apple clang_ compiler on macos Sonoma. It should also work with other _clang_ distributions and with _gcc_. However, _clang_ will produce faster executables as we have not optimized our code for _gcc_.

_CoBarS_ uses several C++ 20 features, so make sure to use a compatible C++ implementation by issueing the compiler option -std=c++20 (or higher).

Parallelization is facilitated by `std::thread` from the _C++ Standard Library_. So you have to use the compiler option `-pthread`. 

Optimization flags like `-O3` or even `-Ofast` are certainly a good idea. I also found that using `-flto` can make a measurable difference (but it ramps up compile time).

With _clang_ as compiler you also have to issue `-fenable-matrix` to enable the clang matrix extension.

See also the example programs in the directories `Example_RandomClosedPolygon`, `Example_Sample_Binned`, and `Example_ConfidenceSample` for usage examples and more detailed compilation instructions.

See also [CoBarSLink](https://github.com/HenrikSchumacher/CoBarSLink) for a more user-friendly _Mathematica_ package.
