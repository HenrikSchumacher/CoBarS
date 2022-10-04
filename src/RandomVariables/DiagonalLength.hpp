#pragma once

#define CLASS DiagonalLength
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
        Real r2 = static_cast<Real>(0);
        
        const SpacePoints_T & p = C.SpaceCoordinates();
        
        const Int last_vertex = C.EdgeCount()/2;
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            const Real delta = p(last_vertex,i) - p(0,i);
            r2 += delta * delta;
        }
        
        return std::sqrt(r2);
    }
    
    virtual Real MinValue( const CyclicSampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const CyclicSampler_T & C ) const override
    {
        const Weights_T & r = C.EdgeLengths();
        
        Real L = static_cast<Real>(0);
        
        const Int last_vertex = C.EdgeCount()/2;
        
        for( Int k = 0; k < last_vertex; ++k )
        {
            L += r[k];
        }
        
        return L;
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return TO_STD_STRING(CLASS);
    }
};

#undef BASE
#undef CLASS
