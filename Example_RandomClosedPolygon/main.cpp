
#include <iostream>
#include "CycleSampler.hpp"

using namespace Tools;
using namespace Tensors;

int main(int argc, const char * argv[])
{
    using Real = float64_t;
    using Int  =   int64_t;

    constexpr Int d            = 3; // Dimensions of the ambient space has to be a compile-time constant.
    const     Int edge_count   = 8;
    const     Int sample_count = 10000000;
    const     Int thread_count = 8; // 0 means "automatic"
    
    // Create an instance of the cycle sampler.
    CycleSampler::Sampler<d,Real,Int>             S (edge_count);
    CycleSampler::MomentPolytopeSampler<Real,Int> M (edge_count);

    // Create containers for the data samples.
    Tensor3<Real,Int> x      ( sample_count, d, edge_count ); // unit edge vectors of open polygons
    Tensor2<Real,Int> w      ( sample_count, d             ); // conformal barycenters
    Tensor3<Real,Int> y      ( sample_count, d, edge_count ); // unit edge vectors of closed polygons
    Tensor1<Real,Int> K      ( sample_count                ); // sample weights for the Pol space
    Tensor1<Real,Int> K_quot ( sample_count                ); // sample weights for the quotient space

    print("");
    print("Settings:");
    
    S.Settings().PrintStats();

    print("");
    valprint("sample_count",sample_count);
    valprint("thread_count",thread_count);
    print("");

    tic("CycleSampler::Sampler::RandomClosedPolygons");
        S.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc("CycleSampler::Sampler::RandomClosedPolygons");
    
    Tensor3<Real,Int> p ( sample_count, d, edge_count+1 ); // vertex positions of polygon
    
    tic("CycleSampler::MomentPolytopeSampler::RandomClosedPolygons");
        M.RandomClosedPolygons(
            p.data(), sample_count, thread_count
        );
    toc("CycleSampler::MomentPolytopeSampler::RandomClosedPolygons");

    
    
    print("");

    valprint("last K     ", K     [sample_count-1], 16 );
    valprint("last K_quot", K_quot[sample_count-1], 16 );

    print("");

    return 0;
}

