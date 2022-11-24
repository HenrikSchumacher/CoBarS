#pragma once

#define CLASS Gyradius
#define BASE  RandomVariable<AmbDim,Real,Int>

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public BASE
{
public:
    
    using Sampler_T         = typename BASE::Sampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    using Weights_T         = typename BASE::Weights_T;
    
    CLASS() = default;
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

protected:
    
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        Real r2 = static_cast<Real>(0);
        
        const Int n             = C.EdgeCount();
        const SpacePoints_T & p = C.SpaceCoordinates();
        
        for( Int k = 0; k < n; ++k )
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                r2 += p(k,i) * p(k,i);
            }
        }
        
        return std::sqrt( r2/n );
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        return Total(C.EdgeLengths())/std::sqrt(C.EdgeCount());
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return TO_STD_STRING(CLASS);
    }
};
        
#undef BASE
#undef CLASS
