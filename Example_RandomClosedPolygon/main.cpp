
#include <iostream>
#include "CycleSampler.hpp"

using namespace Tools;
using namespace Tensors;
using namespace CycleSampler;

int main(int argc, const char * argv[])
{
    
    using Real = double;
    using Int  = int32_t;
    
    constexpr Int d            = 3;
    const     Int edge_count   = 8;
    const     Int sample_count = 10000000;
    const     Int thread_count = 8; // 0 means "automatic"

    
    using Sampler_T = Sampler <d,Real,Int>;
    
    Sampler<d,Real,Int> C (edge_count);

    Tensor3<Real,Int> x      ( sample_count, d, edge_count, 0. );
    Tensor2<Real,Int> w      ( sample_count, d            , 0. );
    Tensor3<Real,Int> y      ( sample_count, d, edge_count, 0. );
    Tensor1<Real,Int> K      ( sample_count               , 0. );
    Tensor1<Real,Int> K_quot ( sample_count               , 0. );

    print("");
    print("Settings:");
    C.Settings().PrintStats();

    print("");
    valprint("sample_count",sample_count);
    valprint("thread_count",thread_count);
    valprint("edge_count  ",edge_count  );
    print("");

    tic("RandomClosedPolygons");
        C.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc("RandomClosedPolygons");

    print("");

    valprint( "last K_out", K[sample_count-1], 16 );

    print("");

    return 0;
}

