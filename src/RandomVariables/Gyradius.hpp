#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class Gyradius : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Sampler_T = Sampler<AmbDim,Real,Int>;
    
    Gyradius() = default;
    
    virtual ~Gyradius() override = default;
    
    __ADD_CLONE_CODE__(Gyradius)

protected:
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        Real r2 = static_cast<Real>(0);
        
        const Int n                   = C.EdgeCount();
        const Real * restrict const p = C.SpaceCoordinates();
        
        for( Int k = 0; k < n; ++k )
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                r2 += p[AmbDim*k+i] * p[AmbDim*k+i];
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
        const Int n = C.EdgeCount();
        
        const Real * restrict r = C.EdgeLengths();
        
        Real sum = 0;
        
        for( Int i = 0; i < n; ++i )
        {
            sum += r[i];
        }
        
        return sum / std::sqrt(C.EdgeCount());
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return "Gyradius";
    }
};
