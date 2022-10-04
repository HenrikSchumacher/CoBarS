#pragma once

#define CLASS GyradiusP
#define BASE  RandomVariable<AmbDim,Real,Int>

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public BASE
{
public:
    
    using CyclicSampler_T   = typename BASE::CyclicSampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    using Weights_T         = typename BASE::Weights_T;
    
    CLASS( const Real exponent_ )
    :   exponent( exponent_ )
    {}
    
    // Copy constructor
    CLASS( const CLASS & other )
    :   exponent(other.exponent)
    {}

    // Move constructor
    CLASS( CLASS && other ) noexcept
    :
        exponent(other.exponent)
    {}
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

protected:
    
    const Real exponent = 2;
    
    virtual Real operator()( const CyclicSampler_T & C ) const override
    {
        Real sum = static_cast<Real>(0);
        Real r2  = static_cast<Real>(0);
        
        const Real power = exponent/2;
        
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
                
                sum+= std::pow(r2,power);
            }
        }
        
        return std::pow( sum/(edge_count * edge_count), static_cast<Real>(1)/exponent );
    }
    
    virtual Real MinValue( const CyclicSampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const CyclicSampler_T & C ) const override
    {
        return Total(C.Omega())/std::pow(C.EdgeCount(),static_cast<Real>(1)/exponent);
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return TO_STD_STRING(CLASS)+"("+ToString(exponent)+")";
    }
    
    virtual std::string ClassName() const override
    {
        return TO_STD_STRING(CLASS)+"<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
    }
};
    
#undef BASE
#undef CLASS
