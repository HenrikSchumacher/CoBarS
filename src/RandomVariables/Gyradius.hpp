#pragma once

namespace CycleSampler
{

#define CLASS Gyradius
    
    template<int AmbDim, typename Real = double, typename Int = long long>
    class CLASS : public RandomVariable<AmbDim,Real,Int>
    {
    private:
        
        using Base_T            = RandomVariable<AmbDim,Real,Int>;
        
    public:
        
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
            
            const Int n             = C.EdgeCount();
            const SpacePoints_T & p = C.SpaceCoordinates();
            
            for( Int k = 0; k < n; ++k )
            {
                for( Int i = 0; i < AmbDim; ++i )
                {
                    r2 += p[k][i] * p[k][i];
                }
            }
            
            return std::sqrt( r2/n );
        }
        
        virtual Real MinValue( const Sampler_T & C ) const override
        {
            return static_cast<Real>(0);
        }
        
        virtual Real MaxValue( const Sampler_T & C ) const override
        {
            return Total(C.EdgeLengths())/std::sqrt(C.EdgeCount());
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
    
}
