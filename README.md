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
    
    constexpr int d = 3; // Dimensions of the ambient space has to be a compile-time constant.
    const int edge_count = 64;
            
    // Create an instance of the cycle sampler.
    // This assumes that equilateral polygons shall be created.
    CoBarS::Sampler<d,double,int> S (edge_count);
    
Then generate the desired edge weights and preallocate memory for the outputs
    
    
    const int sample_count = 10000000;
    
    // The edge vectors of open polygons,
    // to be interpreted as 3-tensor of size `sample_count` x `edge_count` x `d`.
    std::vector<double> x ( sample_count * edge_count * d );
    
    // The conformal barycenters,
    // to be interpreted as 2-tensor of size `sample_count` x `d`.
    std::vector<double> w ( sample_count * d );
    
    // The unit edge vectors of closed polygons,
    // to be interpreted as 3-tensor of `size sample_count` x `edge_count` x `d`.
    std::vector<double> y ( sample_count * edge_count * d );
    
    // The sample weights for the Pol space.
    std::vector<double> K ( sample_count );
    
    // The sample weights for the quotient space.
    std::vector<double> K_quot ( sample_count );
    
    const int thread_count = 8;
    
    // Do the actual sampling.
    S.RandomClosedPolygons(
         &x[0], &w[0], &y[0], &K[0], &K_quot[0], sample_count, thread_count
    );


See also the example programs in the directories Example_RandomClosedPolygon and Example_Sample_Binned for usage examples and compilation instructions.

See also [CoBarSLink](https://github.com/HenrikSchumacher/CoBarSLink) for a more user-friendly _Mathematica_ package.
