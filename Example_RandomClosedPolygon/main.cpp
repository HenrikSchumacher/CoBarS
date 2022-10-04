
#include <iostream>
#include "CyclicSampler.hpp"

using namespace Tools;
using namespace Tensors;

int main(int argc, const char * argv[])
{

    using Real = float64_t;
    using Int  = int32_t;
    
    // insert code here...
    print("Hello, World!");
    constexpr Int d            = 3;
    const     Int edge_count   = 8;
    const     Int sample_count = 10000000;
    const     Int thread_count = 8; // 0 means "automatic"
    
    CyclicSampler::CyclicSampler<d,Real,Int> C (edge_count);

    Tensor3<Real,Int> x      ( sample_count, d, edge_count, 0. );
    Tensor2<Real,Int> w      ( sample_count, d            , 0. );
    Tensor3<Real,Int> y      ( sample_count, d, edge_count, 0. );
    Tensor1<Real,Int> K      ( sample_count               , 0. );
    Tensor1<Real,Int> K_quot ( sample_count               , 0. );

//    print( "Omega = " + C->Omega().ToString() );
//    print( "Rho   = " + C->Rho().ToString() );

    print("");
    print("Settings:");
    C.Settings().PrintStats();

    print("");

    tic("RandomClosedPolygons");
        C.RandomClosedPolygons(
            x.data(),
            w.data(),
            y.data(),
            K.data(),
            K_quot.data(),
            sample_count,
            thread_count
        );
    toc("RandomClosedPolygons");

    print("");

    valprint( "last K_out", K[sample_count-1], 16 );

//    tic("SampleHOMFLY");
//    std::map<std::string, std::tuple<double,double,double>> data = C->SampleHOMFLY( 160, thread_count );
//    toc("SampleHOMFLY");

    print("");

    return 0;
}

