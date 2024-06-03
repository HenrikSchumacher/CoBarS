#include <iostream>
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
    print("Hello, this small test program compares the runtimes of CoBarS::Sampler with various settings to the Action Angle Method (AMM) and the Progressive Action Angle Method (PAMM). Moreover, I use it to detect compile errors in all code paths.");
    
    using Real = double;
    using Int  = int_fast32_t;

    constexpr Int d            = 3; // Dimensions of the ambient space has to be a compile-time constant.
    const     Int edge_count   = 28;
//    const     Int edge_count   = 64;
    const     Int sample_count = 10000000;
    const     Int thread_count = 8;

    
    // Create an instance of the cycle sampler.
    CoBarS::Sampler<d,Real,Int,MT64,false,false>            S_MT_0          (edge_count);
    CoBarS::Sampler<d,Real,Int,MT64,false,true >            S_MT_1          (edge_count);
    CoBarS::Sampler<d,Real,Int,MT64,true ,false>            S_MT_vec_0      (edge_count);
    CoBarS::Sampler<d,Real,Int,MT64,true ,true >            S_MT_vec_1      (edge_count);
    
    CoBarS::Sampler<d,Real,Int,PCG64,false,false>           S_PCG64_0       (edge_count);
    CoBarS::Sampler<d,Real,Int,PCG64,false,true >           S_PCG64_1       (edge_count);
    CoBarS::Sampler<d,Real,Int,PCG64,true ,false>           S_PCG64_vec_0   (edge_count);
    CoBarS::Sampler<d,Real,Int,PCG64,true ,true >           S_PCG64_vec_1   (edge_count);
    
    CoBarS::Sampler<d,Real,Int,Xoshiro256Plus,false,false>  S_Xoshiro_0     (edge_count);
    CoBarS::Sampler<d,Real,Int,Xoshiro256Plus,false,true >  S_Xoshiro_1     (edge_count);
    CoBarS::Sampler<d,Real,Int,Xoshiro256Plus,true ,false>  S_Xoshiro_vec_0 (edge_count);
    CoBarS::Sampler<d,Real,Int,Xoshiro256Plus,true ,true >  S_Xoshiro_vec_1 (edge_count);
    
    CoBarS::Sampler<d,Real,Int,WyRand,false,false>  S_WyRand_0     (edge_count);
    CoBarS::Sampler<d,Real,Int,WyRand,false,true >  S_WyRand_1     (edge_count);
    CoBarS::Sampler<d,Real,Int,WyRand,true ,false>  S_WyRand_vec_0 (edge_count);
    CoBarS::Sampler<d,Real,Int,WyRand,true ,true >  S_WyRand_vec_1 (edge_count);

    
    // The boolean in the AAM::Sampler template stands for progressive (PAAM) or not (AAM).
    
    AAM::Sampler<Real,Int,MT64,false>           M_MT64_0    (edge_count);
    AAM::Sampler<Real,Int,MT64,true >           M_MT64_1    (edge_count);
    
    AAM::Sampler<Real,Int,PCG64,false>          M_PCG64_0   (edge_count);
    AAM::Sampler<Real,Int,PCG64,true >          M_PCG64_1   (edge_count);
    
    AAM::Sampler<Real,Int,Xoshiro256Plus,false> M_Xoshiro_0 (edge_count);
    AAM::Sampler<Real,Int,Xoshiro256Plus,true > M_Xoshiro_1 (edge_count);
    
    AAM::Sampler<Real,Int,WyRand,false> M_WyRand_0 (edge_count);
    AAM::Sampler<Real,Int,WyRand,true > M_WyRand_1 (edge_count);

    // Create containers for the data samples.
    Tensor3<Real,Int> x      ( sample_count, edge_count, d ); // unit edge vectors of open polygons
    Tensor2<Real,Int> w      ( sample_count, d             ); // conformal barycenters
    Tensor3<Real,Int> y      ( sample_count, edge_count, d ); // unit edge vectors of closed polygons
    Tensor1<Real,Int> K      ( sample_count                ); // sample weights for the Pol space
    Tensor1<Real,Int> K_quot ( sample_count                ); // sample weights for the quotient space

    print("");
    print("Settings:");
    
    S_Xoshiro_vec_0.Settings().PrintStats();

    print("");
    valprint("edge_count  ",edge_count  );
    valprint("sample_count",sample_count);
    valprint("thread_count",thread_count);
    print("");
    
    
    
    tic(S_MT_0.ClassName());
        S_MT_0.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_MT_0.ClassName());

    tic(S_MT_1.ClassName());
        S_MT_1.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_MT_1.ClassName());

    tic(S_MT_vec_0.ClassName());
        S_MT_vec_0.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_MT_vec_0.ClassName());

    tic(S_MT_vec_1.ClassName());
        S_MT_vec_1.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_MT_vec_1.ClassName());


    tic(S_PCG64_0.ClassName());
        S_PCG64_0.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_PCG64_0.ClassName());

    tic(S_PCG64_1.ClassName());
        S_PCG64_1.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_PCG64_1.ClassName());

    tic(S_PCG64_vec_0.ClassName());
        S_PCG64_vec_0.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_PCG64_vec_0.ClassName());

    tic(S_PCG64_vec_1.ClassName());
        S_PCG64_vec_1.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_PCG64_vec_1.ClassName());

    
    tic(S_WyRand_0.ClassName());
        S_WyRand_0.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_WyRand_0.ClassName());

    tic(S_WyRand_1.ClassName());
        S_WyRand_1.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_WyRand_1.ClassName());

    tic(S_WyRand_vec_0.ClassName());
        S_WyRand_vec_0.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_WyRand_vec_0.ClassName());

    tic(S_WyRand_vec_1.ClassName());
        S_WyRand_vec_1.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_WyRand_vec_1.ClassName());
    
    
    tic(S_Xoshiro_0.ClassName());
        S_Xoshiro_0.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_Xoshiro_0.ClassName());

    tic(S_Xoshiro_1.ClassName());
        S_Xoshiro_1.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_Xoshiro_1.ClassName());

    tic(S_Xoshiro_vec_0.ClassName());
        S_Xoshiro_vec_0.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_Xoshiro_vec_0.ClassName());

    tic(S_Xoshiro_vec_1.ClassName());
        S_Xoshiro_vec_1.CreateRandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_Xoshiro_vec_1.ClassName());
    


    
    print("");
    
    Tensor3<Real,Int> p ( sample_count, edge_count+1, d ); // vertex positions of polygon

    tic(M_MT64_0.ClassName());
        M_MT64_0.CreateRandomClosedPolygons( p.data(), sample_count, thread_count );
    toc(M_MT64_0.ClassName());

    tic(M_MT64_1.ClassName());
        M_MT64_1.CreateRandomClosedPolygons( p.data(), sample_count, thread_count );
    toc(M_MT64_1.ClassName());
    
    
    tic(M_PCG64_0.ClassName());
        M_PCG64_0.CreateRandomClosedPolygons( p.data(), sample_count, thread_count );
    toc(M_PCG64_0.ClassName());
    
    tic(M_PCG64_1.ClassName());
        M_PCG64_1.CreateRandomClosedPolygons( p.data(), sample_count, thread_count );
    toc(M_PCG64_1.ClassName());
    
    
    tic(M_WyRand_0.ClassName());
        M_WyRand_0.CreateRandomClosedPolygons( p.data(), sample_count, thread_count );
    toc(M_WyRand_0.ClassName());
    
    tic(M_WyRand_1.ClassName());
        M_WyRand_1.CreateRandomClosedPolygons( p.data(), sample_count, thread_count );
    toc(M_WyRand_1.ClassName());
    
    
    tic(M_Xoshiro_0.ClassName());
        M_Xoshiro_0.CreateRandomClosedPolygons( p.data(), sample_count, thread_count );
    toc(M_Xoshiro_0.ClassName());
    
    tic(M_Xoshiro_1.ClassName());
        M_Xoshiro_1.CreateRandomClosedPolygons( p.data(), sample_count, thread_count );
    toc(M_Xoshiro_1.ClassName());

    
    
    print("");

    valprint("last K     ", K     [sample_count-1], 16 );
    valprint("last K_quot", K_quot[sample_count-1], 16 );

    print("");


    return 0;
}
