#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class EdgeQuotientSpaceSamplingWeight : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Sampler_T = Sampler<AmbDim,Real,Int>;
    
    EdgeQuotientSpaceSamplingWeight() = default;
    
    virtual ~EdgeQuotientSpaceSamplingWeight() override = default;
    
    __ADD_CLONE_CODE__(EdgeQuotientSpaceSamplingWeight)

protected:
    
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
