#pragma once

#define CLASS ShiftNorm
#define BASE  RandomVariable<AmbDim,Real,Int>

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public BASE
{
public:
    
    using CyclicSampler_T   = typename BASE::CyclicSampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    
    CLASS() = default;
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

protected:
    
    virtual Real operator()( const CyclicSampler_T & C ) const override
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
    
    virtual Real MinValue( const CyclicSampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const CyclicSampler_T & C ) const override
    {
        return static_cast<Real>(1);
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return TO_STD_STRING(CLASS);
    }
    
    virtual std::string ClassName() const override
    {
        return TO_STD_STRING(CLASS)+"<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
    }
};
        
#undef BASE
#undef CLASS
