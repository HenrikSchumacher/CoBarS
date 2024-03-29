#include <iostream>
#include "CoBarS.hpp"

using namespace Tools;
using namespace Tensors;
using namespace CoBarS;

int main(int argc, const char * argv[])
{
    // Some type aliases to make out lives a bit easier.
    using Real = double;
    using Int  = int_fast32_t;
    
    constexpr Int d            = 3; // Dimensions of the ambient space has to be a compile-time constant.
    const     Int edge_count   = 4;
    const     Int sample_count = 10000000;
    const     Int thread_count = 8;

    // Everything is templated on (i) the dimension of the ambient space, (ii) the floating point type, and (iii) the integer type used, e.g., for indexing.
    
    using SamplerBase_T    = SamplerBase<d,Real,Int>;
    using RandomVariable_T = typename SamplerBase_T::RandomVariable_T;

    print("Test program for routine Sample");
    
    SamplerSettings<Real,Int> opts;
    
    Tensor1<Real,Int> r   ( edge_count, Scalar::One<Real> );
    Tensor1<Real,Int> rho ( edge_count, Scalar::One<Real> );
    Sampler<d,Real,Int,Xoshiro256Plus> S ( r.data(), rho.data(), edge_count, opts );
    
    dump(S.EdgeLengths());
    
    // A list of random variables to sample. We start with an empty list.

//    std::vector< std::shared_ptr<RandomVariableBase_T> > F_list;
    std::vector< std::shared_ptr<RandomVariable_T> > F_list;
    
    // Push as many descendants of RandomVariable_T onto F_list as you like.
    // The nature of runtime polymorphism has it that we have to use smart pointers here...
    F_list.push_back( std::make_shared<ChordLength            <SamplerBase_T>>(0,2));
    F_list.push_back( std::make_shared<ShiftNorm              <SamplerBase_T>>()   );
    F_list.push_back( std::make_shared<EdgeQuotientSpaceSamplingWeight<SamplerBase_T>>() );
    F_list.push_back( std::make_shared<IterationCount         <SamplerBase_T>>()   );
    F_list.push_back( std::make_shared<Gyradius               <SamplerBase_T>>()   );
    F_list.push_back( std::make_shared<HydrodynamicRadius     <SamplerBase_T>>()   );
    F_list.push_back( std::make_shared<EdgeSpaceSamplingWeight<SamplerBase_T>>()   );
    F_list.push_back( std::make_shared<BendingEnergy          <SamplerBase_T>>(2)  );
    F_list.push_back( std::make_shared<MaxAngle               <SamplerBase_T>>()   );
    F_list.push_back( std::make_shared<TotalCurvature         <SamplerBase_T>>()   );

    const Int fun_count    = static_cast<Int>(F_list.size());
    const Int bin_count    = 40;
    const Int moment_count = 3;
    
    // Now we prepare arrays to store the results.
    // S.Sample_Binned will _add into_ these arrays, hence we have to make sure that they are initialized appropriately by zeroes.
    
    
    // bins is a 3D-array of size 3 x fun_count x bin_count. Entry bins(i,j,k) will store the sampled weighted sum in bin k of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i:
    // i == 0: Naive weighting. This is actually the wrong weighting.
    // i == 1: Appropriate reweighting for the total space (without modding out SO(3)).
    // i == 2: Appropriate reweighting for the quotient space (after modding out SO(3)).
    Tensor3<Real,Int> bins    ( 3, fun_count, bin_count,    Real(0) );
    
    // moments is a 3D-array of size 3 x fun_count x bin_count. Entry moments(i,j,k) will store the sampled weighted k-th moment of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i (see above).
    
    Tensor3<Real,Int> moments ( 3, fun_count, moment_count, Real(0) );
    
    // Specify the range for binning: For j-th function in F_list, the range from ranges(j,0) to ranges(j,1) will be divided into bin_count bins.
    Tensor2<Real,Int> ranges  ( fun_count, 2 );

    // The user is supposed to provide meaningful ranges. Some rough guess might be obtained by calling the random variables on the prepared Sampler_T C.
    
    for( Int j = 0; j < fun_count; ++j )
    {
        ranges(j,0) = F_list[j]->MinValue(S);
        ranges(j,1) = F_list[j]->MaxValue(S);
    }

    print("");
    
    // Print the settings.
    print("Settings:");
    S.Settings().PrintStats();

    print("");

    // Perform the actual sampling.
    // The interface operates via raw pointers for more flexibility.
    tic("BinnedSample");
        S.BinnedSample(
            bins.data(),
            bin_count,
            moments.data(),
            moment_count,
            ranges.data(),
            F_list,
            sample_count,
            thread_count
        );
    toc("BinnedSample");
    
    // C.Sample adds into the output arrays, but it does _NOT_ normalize bins and moments. This way we can add further results into them later; we can also simply use them in a further call to C.Sample_Binned.
    
    // Get normalized bins.
    S.NormalizeBinnedSamples(
       bins.data(),    bin_count,
       moments.data(), moment_count,
       fun_count
    );
    
    // Plot a very simplistic histogram
    print("");
    print("Histogram for variable "+F_list[0]->Tag()+" (naive weighting):");
    print("");
    print( "+ <--- " + ToString( ranges(0,0) ) );
    for( Int i = 0; i < bin_count; ++i )
    {
        print( "|"+std::string( static_cast<Int>(bins(0,0,i)*500), '#') );
    }
    print( "+ <--- " + ToString( ranges(0,1) ) );
    print("");
    
    // Plot a very simplistic histogram
    print("");
    print("Histogram for variable "+F_list[0]->Tag()+" (Pol space weighting):");
    print("");
    print( "+ <--- " + ToString( ranges(0,0) ) );
    for( Int i = 0; i < bin_count; ++i )
    {
        print( "|"+std::string( static_cast<Int>(bins(1,0,i)*500), '#') );
    }
    print( "+ <--- " + ToString( ranges(0,1) ) );
    print("");
    
    // Plot a very simplistic histogram
    print("");
    print("Histogram for variable "+F_list[0]->Tag()+" (Quotient space weighting):");
    print("");
    print( "+ <--- " + ToString( ranges(0,0) ) );
    for( Int i = 0; i < bin_count; ++i )
    {
        print( "|"+std::string( static_cast<Int>(bins(2,0,i)*500), '#') );
    }
    print( "+ <--- " + ToString( ranges(0,1) ) );
    print("");

    
    return 0;
}
