#pragma once

#define CLASS HydrodynamicRadius

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

    static constexpr Real eps = std::numeric_limits<Real>::min();
    
protected:
    
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        Real sum = static_cast<Real>(0);
        Real r2  = static_cast<Real>(0);
        
        const Int n             = C.EdgeCount();
        const SpacePoints_T & p = C.SpaceCoordinates();
        
        for( Int k = 0; k < n; ++k )
        {
            for( Int l = k+1; l < n; ++l )
            {
                r2 = static_cast<Real>(0);
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real delta = p(l,i) - p(k,i);
                    
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
        return C.EdgeLengths().Total();
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return TO_STD_STRING(CLASS);
    }
};

#undef BASE
#undef CLASS
