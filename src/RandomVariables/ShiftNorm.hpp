#pragma once

#define CLASS ShiftNorm
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
        return TO_STD_STRING(CLASS);
    }
};
        
#undef BASE
#undef CLASS
