#pragma once

#define CLASS DiagonalLength

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Base_T            = RandomVariable<AmbDim,Real,Int>;
    using Sampler_T         = typename Base_T::Sampler_T;
    using SpherePoints_T    = typename Base_T::SpherePoints_T;
    using SpacePoints_T     = typename Base_T::SpacePoints_T;
    using Weights_T         = typename Base_T::Weights_T;
    
    CLASS() = default;
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

protected:
    
    virtual Real operator()( const Sampler_T & C ) const override
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
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
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
