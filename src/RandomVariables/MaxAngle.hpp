#pragma once

#define CLASS MaxAngle
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
        
        const Int n              = C.EdgeCount();
        const SpherePoints_T & y = C.EdgeCoordinates();

        Real max_angle = static_cast<Real>(0);
        
        {
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( y[n-1], y[0] );
            
            max_angle = std::max(max_angle,phi);
        }
        
        for( Int k = 0; k < n-1; ++k )
        {
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( y[k], y[k+1] );
            
            max_angle = std::max(max_angle,phi);
        }
        
        return max_angle;
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(M_PI);
    }
    
public:
    
    virtual std::string Tag() const override
    {
        return TO_STD_STRING(CLASS);
    }
};
    
#undef BASE
#undef CLASS
