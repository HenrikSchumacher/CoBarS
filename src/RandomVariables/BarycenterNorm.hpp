#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class BarycenterNorm : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Sampler_T = Sampler<AmbDim,Real,Int>;
    
    BarycenterNorm() = default;
    
    virtual ~BarycenterNorm() override = default;
    
    __ADD_CLONE_CODE__(BarycenterNorm)

protected:
    
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        const Real * restrict const p = C.SpaceCoordinates();
        
        const Real * restrict const r = C.EdgeLengths();
        
        const Int n = C.EdgeCount();
        
        Real b[AmbDim];
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            b[i] = r[n-1] * ( p[AmbDim*n+i] + p[AmbDim*0+i] );
        }
        
        for( Int k = 0; k < n-1; ++ k )
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                b[i] += r[i] * ( p[AmbDim*k+i] + p[AmbDim*(k+1)+i] );
            }
        }
        
        Real r2 = static_cast<Real>(0);
        
        Real sum = 0;
        
        for( Int k = 0; k < n; ++k )
        {
            sum += r[k];
        }
        
        const Real factor = static_cast<Real>(0.5)/sum;
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            b[i] *= factor;
            r2 += b[i] * b[i];
        }
        
        return std::sqrt( r2 );
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        return Total(C.Omega());
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return "BarycenterNorm";
    }
};
