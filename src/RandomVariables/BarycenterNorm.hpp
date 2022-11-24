#pragma once

#define CLASS BarycenterNorm
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
        const SpacePoints_T & p = C.SpaceCoordinates();
        
        const Weights_T & r = C.EdgeLengths();
        
        const Int edge_count = C.EdgeCount();
        
        Real b[AmbDim];
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            b[i] = r(edge_count-1) * (p(edge_count,i)+p(0,i));
        }
        
        for( Int k = 0; k < edge_count-1; ++ k )
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                b[i] += r(i) * ( p(k,i) + p(k+1,i) );
            }
        }
        
        Real r2 = static_cast<Real>(0);
        
        const Real factor = static_cast<Real>(0.5)/r.Total();
        
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
        return TO_STD_STRING(CLASS);
    }
};
    
#undef BASE
#undef CLASS
