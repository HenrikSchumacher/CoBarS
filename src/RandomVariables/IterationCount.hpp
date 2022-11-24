#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class IterationCount : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Sampler_T = Sampler<AmbDim,Real,Int>;
    
    IterationCount() = default;
    
    virtual ~IterationCount() override = default;
    
    __ADD_CLONE_CODE__(IterationCount)

protected:
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        return C.IterationCount();
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        return C.MaxIterationCount();
    }
    
public:

    virtual std::string Tag() const  override
    {
        return "IterationCount";
    }
};
