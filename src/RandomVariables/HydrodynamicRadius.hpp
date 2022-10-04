#pragma once

#define CLASS HydrodynamicRadius
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

    static constexpr Real eps = std::numeric_limits<Real>::min();
    
protected:
    
    
    virtual Real operator()( const CyclicSampler_T & C ) const override
    {
        Real sum = static_cast<Real>(0);
        Real r2  = static_cast<Real>(0);
        
        const Int edge_count    = C.EdgeCount();
        const SpacePoints_T & p = C.SpaceCoordinates();
        
        for( Int k = 0; k < edge_count; ++k )
        {
            for( Int l = k+1; l < edge_count; ++l )
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
        
        return (edge_count * edge_count)/sum;
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
    
    virtual bool RequiresSpaceCurve() const override
    {
        return true;
    };
    
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
