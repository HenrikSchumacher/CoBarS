#include "CoBarS.hpp"

// This demonstrate how to create an instance of `CoBarS::Sampler` and how to call its method `RandomClosedPolygons` to generate random polygons along with their weights.


int main()
{
    
    // Dimensions of the ambient space has to be a compile-time constant.
    // Must be a compile-time constant
    constexpr int d = 4;
    
    using Real = double;
    using Int  = std::size_t;
    
    const Int edge_count = 64;
    
    using SamplerBase_T    = CoBarS::SamplerBase<d,Real,Int>;
    using Sampler_T        = CoBarS::Sampler    <d,Real,Int,CoBarS::Xoshiro256Plus>;
    using RandomVariable_T = CoBarS::RandomVariable<SamplerBase_T>;
    
    std::vector<Real> r   ( edge_count, 1 );
    std::vector<Real> rho ( edge_count, 1 );

    Sampler_T S ( &r[0], &rho[0], edge_count );
    
    std::vector< std::shared_ptr<RandomVariable_T> > random_vars;
    
    random_vars.push_back(
        std::make_shared<CoBarS::Gyradius<SamplerBase_T>>()
    );
    
    random_vars.push_back(
        std::make_shared<CoBarS::SquaredGyradius<SamplerBase_T>>()
    );
    
    random_vars.push_back(
        std::make_shared<CoBarS::ChordLength<SamplerBase_T>>(0,10)
    );
    
    const double relative_confidence_radius = 0.0005;
    
    std::vector<Real> radii            ( random_vars.size(), relative_confidence_radius );
    
    std::vector<Real> means            ( random_vars.size(), 0 );
    std::vector<Real> sample_variances ( random_vars.size(), 0 );
    std::vector<Real> errors           ( random_vars.size(), 0 );
    
    
    const Int  max_sample_count = 100000000;
    const bool quot_space_Q     = true;
    const Int  thread_count     = 8;
    const Real confidence_level = 0.99;
    const Int  chunk_size       = 1000000;
    const bool relativeQ        = true;
    const bool verboseQ         = true;
    
    S.ConfidenceSample(
        random_vars,
        &means[0],
        &sample_variances[0],
        &errors[0],
        &radii[0],
        max_sample_count,
        quot_space_Q,
        thread_count,
        confidence_level,
        chunk_size,
        relativeQ,
        verboseQ
    );
    
    Tools::print("");
    
    dump( means );
    dump( sample_variances );
    dump( errors );
    
    return 0;
}

