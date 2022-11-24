#pragma once

namespace CycleSampler {

#define CLASS RandomVariable
#define BASE  RandomVariableBase<Real,Int>

    template<int AmbDim, typename Real = double, typename Int = long long>
    class CLASS : public BASE
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using SamplerBase_T     = SamplerBase<Real,Int>;
        using Sampler_T         = Sampler<AmbDim,Real,Int>;
        using SpherePoints_T    = typename Sampler_T::SpherePoints_T;
        using SpacePoints_T     = typename Sampler_T::SpacePoints_T;
        using Weights_T         = typename Sampler_T::Weights_T;
        using Vector_T          = typename Sampler_T::Vector_T;
        
        CLASS() = default;
        
        virtual ~CLASS(){}
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            return operator()( dynamic_cast<const Sampler_T &>(C) );
        }
        
        virtual Real operator()( const Sampler_T & C ) const = 0;
        
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return MinValue( dynamic_cast<const Sampler_T &>(C) );
        }
        virtual Real MinValue( const Sampler_T & C ) const = 0;
        
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return MaxValue( dynamic_cast<const Sampler_T &>(C) );
        }
        
        virtual Real MaxValue( const Sampler_T & C ) const = 0;
        
        __ADD_CLONE_CODE_FOR_ABSTRACT_CLASS__(CLASS)

    public:
        
        Int AmbientDimension() const override
        {
            return AmbDim;
        }
        
        virtual std::string Tag() const override = 0;
    };
        
#undef BASE
#undef CLASS
    
} // namespace CycleSampler
