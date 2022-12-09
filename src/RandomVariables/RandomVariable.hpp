#pragma once

namespace CycleSampler
{

    template<int AmbDim, typename Real = double, typename Int = long long>
    class RandomVariable
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using Sampler_T         = Sampler<AmbDim,Real,Int>;
        using SpherePoints_T    = typename Sampler_T::SpherePoints_T;
        using SpacePoints_T     = typename Sampler_T::SpacePoints_T;
        using Weights_T         = typename Sampler_T::Weights_T;
        using Vector_T          = typename Sampler_T::Vector_T;
        
        RandomVariable() = default;
        
        virtual ~RandomVariable(){}
        
        virtual Real operator()( const Sampler_T & C ) const = 0;
        
        virtual Real MinValue( const Sampler_T & C ) const = 0;
        
        virtual Real MaxValue( const Sampler_T & C ) const = 0;
        
        __ADD_CLONE_CODE_FOR_BASE_CLASS__(RandomVariable)

    public:
        
        Int AmbientDimension() const
        {
            return AmbDim;
        }
        
        virtual std::string Tag() const = 0;
        
    }; // RandomVariable
    
} // namespace CycleSampler
