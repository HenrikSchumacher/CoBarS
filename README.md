# CoBarS - Conformal Barycenter Sampling

by Jason Cantarella and Henrik Schumacher

A header-only C++ library for Monte-Carlo sampling of cylic polygons with prescribed edge lengths.

# Installation

Please clone with

    git clone --recurse-submodules git@github.com:HenrikSchumacher/CoBarS.git

to load also all submodules. If you forgot to do that, you can also run the following afterwards:

    git submodule update --init --recursive
    
# Usage

Just include the header `CoBarS.hpp` into your C++ program via

    #import "CoBarS.hpp"    
    
First create a `CoBarS::Sampler` object. 
    
    constexpr int d            = 3; // Dimensions of the ambient space has to be a compile-time constant.
    const     int edge_count   = 64;
            
    // Create an instance of the cycle sampler.
    CoBarS::Sampler<d,double,int> S(edge_count);
    
Then generate the desired edge weights and preallocate memory for the outputs
    
    
    const     int sample_count = 10000000;
    
    std::vector<double> w
    
    
    const     int thread_count = 8;
    
    
            S_MT_0.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );


See also the examples programs in the directories Example_RandomClosedPolygon and Example_Sample_Binned for usage examples.

See also [CoBarSLink](https://github.com/HenrikSchumacher/CoBarSLink) for a more user-friendly _Mathematica_ package.
