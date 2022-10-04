#pragma once

#define CLASS IterationCount
#define BASE  RandomVariable<AmbDim,Real,Int>

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public BASE
{
public:
    
    using CyclicSampler_T   = typename BASE::CyclicSampler_T;
    
    CLASS() = default;
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

protected:
    
    virtual Real operator()( const CyclicSampler_T & C ) const override
    {
        return C.IterationCount();
    }
    
    virtual Real MinValue( const CyclicSampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const CyclicSampler_T & C ) const override
    {
        return C.MaxIterationCount();
    }
    
public:

    virtual std::string Tag() const  override
    {
        return TO_STD_STRING(CLASS);
    }
};
    
#undef BASE
#undef CLASS
