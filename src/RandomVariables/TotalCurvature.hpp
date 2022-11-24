#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class TotalCurvature : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Sampler_T = Sampler<AmbDim,Real,Int>;
    
    TotalCurvature() = default;
    
    virtual ~TotalCurvature() override = default;
    
    __ADD_CLONE_CODE__(TotalCurvature)

protected:
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        
        const Int n                   = C.EdgeCount();
        const Real * restrict const y = C.EdgeCoordinates();

        Real sum;
        
        {
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( &y[AmbDim*(n-1)], &y[AmbDim*0] );
            
            sum = phi;
        }
        
        for( Int k = 0; k < n-1; ++k )
        {
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( &y[AmbDim*k], &y[AmbDim*(k+1)] );
            
            sum += phi;
        }
        
        return sum;
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        return C.EdgeCount() * static_cast<Real>(M_PI);
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return "TotalCurvature";
    }
};
