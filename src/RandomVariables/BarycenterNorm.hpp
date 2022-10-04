#pragma once

#define CLASS BarycenterNorm
#define BASE  RandomVariable<AmbDim,Real,Int>

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public BASE
{
public:
    
    using CyclicSampler_T   = typename BASE::CyclicSampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    using Weights_T         = typename BASE::Weights_T;
    
    CLASS() = default;
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

protected:
    
    
    virtual Real operator()( const CyclicSampler_T & C ) const override
    {
        const SpacePoints_T & p = C.SpaceCoordinates();
        
        const Int edge_count = C.EdgeCount();
        
        Real b[AmbDim];
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            b[i] = static_cast<Real>(0.5) * (p(0,i) + p(edge_count,i));
        }
        
        for( Int k = 1; k < edge_count; ++ k )
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                b[i] += p(k,i);
            }
        }
        
        Real r2 = static_cast<Real>(0);
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            b[i] /= edge_count;
            r2 += b[i] * b[i];
        }
        
        return std::sqrt( r2 );
    }
    
    virtual Real MinValue( const CyclicSampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const CyclicSampler_T & C ) const override
    {
        return Total(C.Omega());
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
