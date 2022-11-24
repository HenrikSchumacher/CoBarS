#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class ShiftNorm : public RandomVariable<AmbDim,Real,Int>
{
public:
    
<<<<<<< HEAD
    using Sampler_T   = typename BASE::Sampler_T;
=======
    using Sampler_T = Sampler<AmbDim,Real,Int>;
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    
    ShiftNorm() = default;
    
    virtual ~ShiftNorm() override = default;
    
    __ADD_CLONE_CODE__(ShiftNorm)

protected:
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        Real r2 = static_cast<Real>(0);
        
        Real w[AmbDim];
        
        C.WriteShiftVector( &w[0] );
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            r2 += w[i] * w[i];
        }
        
        return std::sqrt( r2 );
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(1);
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return "ShiftNorm";
    }
};
