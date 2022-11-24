#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class MaxAngle : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Sampler_T = Sampler<AmbDim,Real,Int>;
    
    MaxAngle() = default;
    
    virtual ~MaxAngle() override = default;
    
    __ADD_CLONE_CODE__(MaxAngle)

protected:
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        
        const Int n                   = C.EdgeCount();
        const Real * restrict const y = C.EdgeCoordinates();

        Real max_angle = static_cast<Real>(0);
        
        {
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( &y[AmbDim*(n-1)], &y[AmbDim*0] );
            
            max_angle = std::max(max_angle,phi);
        }
        
        for( Int k = 0; k < n-1; ++k )
        {
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( &y[AmbDim*k], &y[AmbDim*(k+1)] );
            
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
        return "MaxAngle";
    }
};
