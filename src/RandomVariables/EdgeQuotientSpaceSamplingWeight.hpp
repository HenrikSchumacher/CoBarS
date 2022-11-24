#pragma once

#define CLASS EdgeQuotientSpaceSamplingWeight
#define BASE  RandomVariable<AmbDim,Real,Int>

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public BASE
{
public:
    
    using Sampler_T         = typename BASE::Sampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    
    CLASS() = default;
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

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
        return TO_STD_STRING(CLASS);
    }
};
        
#undef BASE
#undef CLASS
