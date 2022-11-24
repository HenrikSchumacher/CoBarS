#pragma once

#define CLASS EdgeSpaceSamplingWeight
#define BASE  RandomVariable<AmbDim,Real,Int>

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public BASE
{
public:
    
    using Sampler_T   = typename BASE::Sampler_T;
    
    CLASS() = default;
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

protected:
    
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        return C.EdgeSpaceSamplingWeight();
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
//            return static_cast<Real>(1)/( std::pow( C.EdgeCount(), AmbDim-1) );
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
