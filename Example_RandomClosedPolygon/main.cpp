#include <iostream>
#include "CycleSampler.hpp"

using namespace Tools;
using namespace Tensors;
using namespace CycleSampler;

int main(int argc, const char * argv[])
{
    using Real = double;
    using Int  = int_fast32_t;

    constexpr Int d            = 3; // Dimensions of the ambient space has to be a compile-time constant.
    const     Int edge_count   = 24;
    const     Int sample_count = 10000000;
    const     Int thread_count = 8;
    
    // Create an instance of the cycle sampler.
    Sampler<d,Real,Int,MT64,false,false>            S_MT_0          (edge_count);
    Sampler<d,Real,Int,MT64,false,true >            S_MT_1          (edge_count);
    Sampler<d,Real,Int,MT64,true ,false>            S_MT_vec_0      (edge_count);
    Sampler<d,Real,Int,MT64,true ,true >            S_MT_vec_1      (edge_count);
    
    Sampler<d,Real,Int,PCG64,false,false>           S_PCG64_0       (edge_count);
    Sampler<d,Real,Int,PCG64,false,true >           S_PCG64_1       (edge_count);
    Sampler<d,Real,Int,PCG64,true ,false>           S_PCG64_vec_0   (edge_count);
    Sampler<d,Real,Int,PCG64,true ,true >           S_PCG64_vec_1   (edge_count);
    
    Sampler<d,Real,Int,Xoshiro256Plus,false,false>  S_Xoshiro_0     (edge_count);
    Sampler<d,Real,Int,Xoshiro256Plus,false,true >  S_Xoshiro_1     (edge_count);
    Sampler<d,Real,Int,Xoshiro256Plus,true ,false>  S_Xoshiro_vec_0 (edge_count);
    Sampler<d,Real,Int,Xoshiro256Plus,true ,true >  S_Xoshiro_vec_1 (edge_count);
    
    ActionAngleSampler<Real,Int>                    M               (edge_count);

    // Create containers for the data samples.
    Tensor3<Real,Int> x      ( sample_count, d, edge_count ); // unit edge vectors of open polygons
    Tensor2<Real,Int> w      ( sample_count, d             ); // conformal barycenters
    Tensor3<Real,Int> y      ( sample_count, d, edge_count ); // unit edge vectors of closed polygons
    Tensor1<Real,Int> K      ( sample_count                ); // sample weights for the Pol space
    Tensor1<Real,Int> K_quot ( sample_count                ); // sample weights for the quotient space

    print("");
    print("Settings:");
    
    S_Xoshiro_0.Settings().PrintStats();

    print("");
    valprint("edge_count  ",edge_count  );
    valprint("sample_count",sample_count);
    valprint("thread_count",thread_count);
    print("");
    
    
    
    tic(S_MT_0.ClassName());
        S_MT_0.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_MT_0.ClassName());

    tic(S_MT_1.ClassName());
        S_MT_1.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_MT_1.ClassName());

    tic(S_MT_vec_0.ClassName());
        S_MT_vec_0.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_MT_vec_0.ClassName());

    tic(S_MT_vec_1.ClassName());
        S_MT_vec_1.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_MT_vec_1.ClassName());



    tic(S_PCG64_0.ClassName());
        S_PCG64_0.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_PCG64_0.ClassName());

    tic(S_PCG64_1.ClassName());
        S_PCG64_1.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_PCG64_1.ClassName());

    tic(S_PCG64_vec_0.ClassName());
        S_PCG64_vec_0.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_PCG64_vec_0.ClassName());

    tic(S_PCG64_vec_1.ClassName());
        S_PCG64_vec_1.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_PCG64_vec_1.ClassName());

    tic(S_Xoshiro_0.ClassName());
        S_Xoshiro_0.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_Xoshiro_0.ClassName());

    tic(S_Xoshiro_1.ClassName());
        S_Xoshiro_1.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_Xoshiro_1.ClassName());

    tic(S_Xoshiro_vec_0.ClassName());
        S_Xoshiro_vec_0.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_Xoshiro_vec_0.ClassName());

    tic(S_Xoshiro_vec_1.ClassName());
        S_Xoshiro_vec_1.RandomClosedPolygons(
            x.data(), w.data(), y.data(), K.data(), K_quot.data(), sample_count, thread_count
        );
    toc(S_Xoshiro_vec_1.ClassName());

    Tensor3<Real,Int> p ( sample_count, d, edge_count+1 ); // vertex positions of polygon

    tic(M.ClassName());
        M.RandomClosedPolygons(
            p.data(), sample_count, thread_count
        );
    toc(M.ClassName());

    
    
    print("");

    valprint("last K     ", K     [sample_count-1], 16 );
    valprint("last K_quot", K_quot[sample_count-1], 16 );

    print("");


    return 0;
}
