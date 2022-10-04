//
//  main.cpp
//  Example


#define REAL double
#define INT  int

#include <iostream>
#include "CyclicSampler.hpp"

using namespace Tools;
using namespace Tensors;

int main(int argc, const char * argv[]) {
    // insert code here...
    print("Hello, World!");

    std::unique_ptr<CyclicSampler::CyclicSamplerBase<REAL,INT>> C ;


    const INT d            = 3;
    const INT edge_count   = 8;
    const INT sample_count = 10000000;
    const INT thread_count = 8; // 0 means "automatic"

    C = CyclicSampler::MakeCyclicSampler<REAL,INT>(d,edge_count);

    Tensor3<REAL,INT> x      ( sample_count, d, edge_count, 0. );
    Tensor2<REAL,INT> w      ( sample_count, d            , 0. );
    Tensor3<REAL,INT> y      ( sample_count, d, edge_count, 0. );
    Tensor1<REAL,INT> K      ( sample_count               , 0. );
    Tensor1<REAL,INT> K_quot ( sample_count               , 0. );

//    print( "Omega = " + C->Omega().ToString() );
//    print( "Rho   = " + C->Rho().ToString() );

    print("");
    print("Settings:");
    C->Settings().PrintStats();

    print("");

    tic("RandomClosedPolygons");
        C->RandomClosedPolygons(
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
//    std::map<std::string, std::tuple<REAL,REAL,REAL>> data = C->SampleHOMFLY( 160, thread_count );
//    toc("SampleHOMFLY");

    print("");

    return 0;
}
