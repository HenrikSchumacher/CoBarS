#pragma once
template<int AmbDim, typename Real = double, typename Int = long long>
class HydrodynamicRadius : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Sampler_T = Sampler<AmbDim,Real,Int>;
    
    HydrodynamicRadius() = default;
    
    virtual ~HydrodynamicRadius() override = default;
    
    __ADD_CLONE_CODE__(HydrodynamicRadius)

    static constexpr Real eps = std::numeric_limits<Real>::min();
    
protected:
    
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        Real sum = static_cast<Real>(0);
        Real r2  = static_cast<Real>(0);
        
        const Int n                   = C.EdgeCount();
        const Real * restrict const p = C.SpaceCoordinates();
        
        for( Int k = 0; k < n; ++k )
        {
            for( Int l = k+1; l < n; ++l )
            {
                r2 = static_cast<Real>(0);
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real delta = p[AmbDim*l+i] - p[AmbDim*k+i];
                    
                    r2 += delta * delta;
                }
                
                sum+= static_cast<Real>(1)/(std::sqrt(r2) + eps);
            }
        }
        
        return (n * n)/sum;
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        const Int n = C.EdgeCount();
        
        const Real * restrict r = C.EdgeLengths();
        
        Real sum = 0;
        
        for( Int i = 0; i < n; ++i )
        {
            sum += r[i];
        }
        
        return sum;
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return "HydrodynamicRadius";
    }
};
