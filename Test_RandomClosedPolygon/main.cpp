#include <iostream>

#define TOOLS_AGGRESSIVE_INLINING
#define TOOLS_AGGRESSIVE_UNROLLING

#include "CoBarS.hpp"

using namespace Tools;
using namespace Tensors;

// There are multiple pseudorandom number generators available in CoBarS.
// In order to unify the interface, I had to write wrapper classes for them:

// The Mersenne-Twister from the standard library.
using CoBarS::MT64;

// Permuted congruential generator by Melissa O'Neill.
using CoBarS::PCG64;

// A wrapper for wyrand, written by Alain Espinosa.
using CoBarS::WyRand;

// Xoshiro256+ by David Blackman and Sebastiano Vigna.
using CoBarS::Xoshiro256Plus;

// Typically, CoBarS::Xoshiro256Plus is the fastest number generator for doubles, because it creates only the 53 random bits required for a double---not more.
// CoBarS::WyRand is almost as fast; it may have

// See src/MT64.hpp, src/PCG64.hpp, src/WyRand.hpp, and src/Xoshiro256Plus.hpp for the implementation. This should also tell you how to roll out the pseudorandom number generator of your choice.

int main()
{
    print("Hello, this small test program compares the runtimes of CoBarS::Sampler with various settings to the Action Angle Method (AAM) and the Progressive Action Angle Method (PAAM). Moreover, I use it to detect compilation errors in all code paths.");
    
    using Real = double;
    using Int  = std::size_t;

    // Dimensions of the ambient space has to be a compile-time constant.
    constexpr Int d            = 3;
    const     Int edge_count   = 32;

    const     Int sample_count = 10000000;
    const     Int thread_count = 8;
    
    // Whether we want the sampling weights for the quotient space of polygons 
    // modulo SO(d) (`quot_space_Q = true`) or not (`quot_space_Q = false`).
    const     bool quot_space_Q = true;
    
//     Create an instance of the cycle sampler.
    CoBarS::Sampler<d,Real,Int,MT64,false,false>            S_MT_0          (edge_count);
    CoBarS::Sampler<d,Real,Int,MT64,false,true >            S_MT_1          (edge_count);
    CoBarS::Sampler<d,Real,Int,MT64,true ,false>            S_MT_vec_0      (edge_count);
    CoBarS::Sampler<d,Real,Int,MT64,true ,true >            S_MT_vec_1      (edge_count);
    
    CoBarS::Sampler<d,Real,Int,PCG64,false,false>           S_PCG64_0       (edge_count);
    CoBarS::Sampler<d,Real,Int,PCG64,false,true >           S_PCG64_1       (edge_count);
    CoBarS::Sampler<d,Real,Int,PCG64,true ,false>           S_PCG64_vec_0   (edge_count);
    CoBarS::Sampler<d,Real,Int,PCG64,true ,true >           S_PCG64_vec_1   (edge_count);
    
    CoBarS::Sampler<d,Real,Int,WyRand,false,false>          S_WyRand_0      (edge_count);
    CoBarS::Sampler<d,Real,Int,WyRand,false,true >          S_WyRand_1      (edge_count);
    CoBarS::Sampler<d,Real,Int,WyRand,true ,false>          S_WyRand_vec_0  (edge_count);
    CoBarS::Sampler<d,Real,Int,WyRand,true ,true >          S_WyRand_vec_1  (edge_count);
    
    CoBarS::Sampler<d,Real,Int,Xoshiro256Plus,false,false>  S_Xoshiro_0     (edge_count);
    CoBarS::Sampler<d,Real,Int,Xoshiro256Plus,false,true >  S_Xoshiro_1     (edge_count);
    CoBarS::Sampler<d,Real,Int,Xoshiro256Plus,true ,false>  S_Xoshiro_vec_0 (edge_count);
    CoBarS::Sampler<d,Real,Int,Xoshiro256Plus,true ,true >  S_Xoshiro_vec_1 (edge_count);


    
    // The boolean in the AAM::Sampler template stands for progressive (PAAM) or not (AAM).
    
    AAM::Sampler<Real,Int,MT64,false>           M_MT64_0    (edge_count);
    AAM::Sampler<Real,Int,MT64,true >           M_MT64_1    (edge_count);
    
    AAM::Sampler<Real,Int,PCG64,false>          M_PCG64_0   (edge_count);
    AAM::Sampler<Real,Int,PCG64,true >          M_PCG64_1   (edge_count);

    AAM::Sampler<Real,Int,WyRand,false>         M_WyRand_0  (edge_count);
    AAM::Sampler<Real,Int,WyRand,true >         M_WyRand_1  (edge_count);
    
    AAM::Sampler<Real,Int,Xoshiro256Plus,false> M_Xoshiro_0 (edge_count);
    AAM::Sampler<Real,Int,Xoshiro256Plus,true > M_Xoshiro_1 (edge_count);
    
    
    // Create containers for the data samples.
    Tensor3<Real,Int> p ( sample_count, edge_count + 1, d ); // vertex positions of polygons.
    Tensor1<Real,Int> K ( sample_count                    ); // sampling weights

    print("");
    print("Settings:");
    
    S_Xoshiro_vec_0.Settings().PrintStats();

    print("");
    valprint("edge_count  ",edge_count  );
    valprint("sample_count",sample_count);
    valprint("thread_count",thread_count);
    print("");
    
    auto run_CoBarS = [&p,&K,sample_count,quot_space_Q,thread_count]( auto & S )
    {
        tic(S.ClassName());
        
        S.CreateRandomClosedPolygons(
            p.data(), K.data(), sample_count, quot_space_Q, thread_count
        );
        
        toc(S.ClassName());
    };
    
    auto run_AAM = [&p,sample_count,thread_count]( auto & S )
    {
        tic(S.ClassName());
        
        S.CreateRandomClosedPolygons(
            p.data(), sample_count, thread_count
        );
        
        toc(S.ClassName());
    };
    
    
    run_CoBarS(S_MT_1);
    run_CoBarS(S_MT_0);
    run_CoBarS(S_MT_vec_1);
    run_CoBarS(S_MT_vec_0);
    
    print("");
    
    run_CoBarS(S_PCG64_1);
    run_CoBarS(S_PCG64_0);
    run_CoBarS(S_PCG64_vec_1);
    run_CoBarS(S_PCG64_vec_0);
    
    print("");
    
    run_CoBarS(S_WyRand_1);
    run_CoBarS(S_WyRand_0);
    run_CoBarS(S_WyRand_vec_1);
    run_CoBarS(S_WyRand_vec_0);
    
    print("");

    run_CoBarS(S_Xoshiro_1);
    run_CoBarS(S_Xoshiro_0);
    run_CoBarS(S_Xoshiro_vec_1);
    run_CoBarS(S_Xoshiro_vec_0);

    print("");
    
//    run_AAM(M_MT64_0);
//    run_AAM(M_MT64_1);
//    
//    print("");
//
//    run_AAM(M_PCG64_0);
//    run_AAM(M_PCG64_1);
//    
//    print("");
//
//    run_AAM(M_WyRand_0);
//    run_AAM(M_WyRand_1);
//  
//    print("");

    run_AAM(M_Xoshiro_0);
    run_AAM(M_Xoshiro_1);

    return 0;
}
