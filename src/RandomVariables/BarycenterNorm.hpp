#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class BarycenterNorm : public RandomVariable<AmbDim,Real,Int>
{
public:
    
<<<<<<< HEAD
    using Sampler_T         = typename BASE::Sampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    using Weights_T         = typename BASE::Weights_T;
=======
    using Sampler_T = Sampler<AmbDim,Real,Int>;
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    
    BarycenterNorm() = default;
    
    virtual ~BarycenterNorm() override = default;
    
    __ADD_CLONE_CODE__(BarycenterNorm)

protected:
    
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
<<<<<<< HEAD
        const SpacePoints_T & p = C.SpaceCoordinates();
=======
        const Real * restrict const p = C.SpaceCoordinates();
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
        
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
