#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class EdgeQuotientSpaceSamplingWeight : public RandomVariable<AmbDim,Real,Int>
{
public:
    
<<<<<<< HEAD
    using Sampler_T         = typename BASE::Sampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
=======
    using Sampler_T = Sampler<AmbDim,Real,Int>;
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    
    EdgeQuotientSpaceSamplingWeight() = default;
    
    virtual ~EdgeQuotientSpaceSamplingWeight() override = default;
    
    __ADD_CLONE_CODE__(EdgeQuotientSpaceSamplingWeight)

protected:
    
<<<<<<< HEAD
    
=======
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    virtual Real operator()( const Sampler_T & C ) const override
    {
        return C.EdgeQuotientSpaceSamplingWeight();
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(1)/( C.EdgeCount() );
    }
    
public:
    virtual std::string Tag() const  override
    {
        return "EdgeQuotientSpaceSamplingWeight";
    }
};
