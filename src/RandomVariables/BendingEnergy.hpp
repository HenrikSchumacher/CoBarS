#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class BendingEnergy : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Base_T            = RandomVariable<AmbDim,Real,Int>;
    using Sampler_T         = typename Base_T::Sampler_T;
    using SpherePoints_T    = typename Base_T::SpherePoints_T;
    using SpacePoints_T     = typename Base_T::SpacePoints_T;
    using Weights_T         = typename Base_T::Weights_T;
    
    explicit BendingEnergy( const Real p_ )
    :   p( p_ )
    {}
    
    // Copy constructor
    explicit BendingEnergy( const BendingEnergy & other )
    :   p( other.p )
    {}

    // Move constructor
    explicit BendingEnergy( BendingEnergy && other ) noexcept
    :   p( other.p )
    {}
    
    virtual ~BendingEnergy() override = default;
    
    __ADD_CLONE_CODE__(BendingEnergy)

protected:
    
    const Real p;
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        const Int n              = C.EdgeCount();
        const SpherePoints_T & y = C.EdgeCoordinates();
        const Weights_T & r      = C.EdgeLengths();
    
        Real sum;
        
        {
            const Real len = static_cast<Real>(0.5)*(r[n-1]+r[0]);
            
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( y.data(n-1), y.data(0) );
            
            sum = std::pow( phi / len, p ) * len;
        }
        
        for( Int k = 0; k < n-1; ++k )
        {
            const Real len = static_cast<Real>(0.5)*(r[k]+r[k+1]);
            
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( y.data(k), y.data(k+1) );
            
            sum += std::pow( phi / len, p ) * len;
        }
        
        return sum/p;
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        const Int n         = C.EdgeCount();
        const Weights_T & r = C.EdgeLengths();
    
        Real sum;
        
        {
            const Real len = static_cast<Real>(0.5)*(r[n-1]+r[0]);
            
            const Real phi = static_cast<Real>(M_PI);
            
            sum = std::pow( phi / len, p ) * len;
            
        }
        
        for( Int k = 0; k < n-1; ++k )
        {
            const Real len = static_cast<Real>(0.5)*(r[k]+r[k+1]);
            
            const Real phi = static_cast<Real>(M_PI);
            
            sum += std::pow( phi / len, p ) * len;
        }
        
        return sum/p;
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return "BendingEnergy("+ToString(p)+")";
    }
};
    
