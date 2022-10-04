
#include <iostream>
#include "CyclicSampler.hpp"

using namespace Tools;
using namespace Tensors;

int main(int argc, const char * argv[])
{
    // Some type aliases to make out lives a bit easier.
    using Real = float64_t;
    using Int  = int32_t;
    
    constexpr Int AmbDim = 3;

    // Everything is templated on (i) the dimension of the ambient space, (ii) the floating point type, and (iii) the integer type used, e.g., for indexing.
    // In particular, the ambient dimension has to be known at compile time.
    using CyclicSampler_T  = CyclicSampler::CyclicSampler <AmbDim,Real,Int>;
    using RandomVariable_T = CyclicSampler::RandomVariableBase<Real,Int>;
    
    const     Int edge_count   = 8;
    const     Int sample_count = 1000000;
    const     Int thread_count = 8;
    
    print("Test program for routine CyclicSampler::Sample");

    CyclicSampler_T C (edge_count);
    
    // A list of random variables to sample. We start with an empty list.
    std::vector< std::unique_ptr<RandomVariable_T> > F_list;
    
    // Push as many descendants of RandomVariable_T onto F_list as you like.
    // The nature of runtime polymorphism has it that we have to use smart pointers here...
    F_list.push_back( std::make_unique<CyclicSampler::ShiftNorm<AmbDim,Real,Int>>() );
    F_list.push_back( std::make_unique<CyclicSampler::Gyradius<AmbDim,Real,Int>>() );
    F_list.push_back( std::make_unique<CyclicSampler::ChordLength<AmbDim,Real,Int>>(0,2) );
    F_list.push_back( std::make_unique<CyclicSampler::TotalCurvature<AmbDim,Real,Int>>() );

    const Int fun_count    = static_cast<Int>(F_list.size());
    const Int bin_count    = 40;
    const Int moment_count = 3;
    
    // Now we prepare arrays to store the results.
    // C.Sample_Binned will _add into_ these arrays, hence we have to make sure that they are initialized appropriately by zeroes.
    
    
    // bins is a 3D-array of size 3 x fun_count x bin_count. Entry bins(i,j,k) will store the sampled weighted sum in bin k of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i:
    // i == 0: Naive weighting. This is actually the wrong weighting.
    // i == 1: Appropriate reweighting for the total space (without modding out SO(3)).
    // i == 2: Appropriate reweighting for the quotient space (after modding out SO(3)).
    Tensor3<Real,Int> bins    ( 3, fun_count, bin_count,    Real(0) );
    
    // moments is a 3D-array of size 3 x fun_count x bin_count. Entry moments(i,j,k) will store the sampled weighted k-th moment of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i (see above).
    
    Tensor3<Real,Int> moments ( 3, fun_count, moment_count, Real(0) );
    
    // Specify the range for binning: For j-th function in F_list, the range from ranges(j,0) to ranges(j,1) will be devided into bin_count bins.
    Tensor2<Real,Int> ranges  ( fun_count, 2 );

    // The user is supposed to provide meaningful ranges. Some rough guess might be obtained by calling the random variables on the prepared CyclicSampler_T C.
    
    for( Int j = 0; j < fun_count; ++j )
    {
        ranges(j,0) = F_list[j]->MinValue( C );
        ranges(j,1) = F_list[j]->MaxValue( C );
    }

    print("");
    
    // Print the settings.
    print("Settings:");
    C.Settings().PrintStats();

    print("");

    // Perform the actual sampling.
    // The interface operates via raw pointers for more flexibility.
    tic("Sample_Binned");
        C.Sample_Binned(
            bins.data(),
            bin_count,
            moments.data(),
            moment_count,
            ranges.data(),
            F_list,
            sample_count,
            thread_count
        );
    toc("Sample_Binned");
    
    // C.Sample adds into the output arrays, but it does _NOT_ normalize bins and moments. This way we can add futher results into them later; we can also simply use them in a further call to C.Sample.
    
    // In order to get normalized bins.
    
    // The interface operates via raw pointers for more flexibility.
    C.NormalizeBinnedSamples(
       bins.data(),
       bin_count,
       moments.data(),
       moment_count,
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

    return 0;
}

