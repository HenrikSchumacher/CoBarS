#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class ShiftNorm : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Base_T            = RandomVariable<AmbDim,Real,Int>;
    using Sampler_T         = typename Base_T::Sampler_T;
    
    ShiftNorm() = default;
    
    virtual ~ShiftNorm() override = default;
    
    __ADD_CLONE_CODE__(ShiftNorm)

protected:
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        Real r2 = static_cast<Real>(0);
        
        Real w[AmbDim];
        
        C.WriteShiftVector( &w[0] );
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            r2 += w[i] * w[i];
        }
        
        return std::sqrt( r2 );
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(1);
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return "ShiftNorm";
    }
};
