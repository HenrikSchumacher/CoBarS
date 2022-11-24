#pragma once

#define CLASS EdgeSpaceSamplingWeight

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Base_T      = RandomVariable<AmbDim,Real,Int>;
    using Sampler_T   = typename Base_T::Sampler_T;
    
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
