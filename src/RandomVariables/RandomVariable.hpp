#pragma once

namespace CycleSampler
{

    template<int AmbDim, typename Real = double, typename Int = long long>
    class RandomVariable : public RandomVariableBase<Real,Int>
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
        
        RandomVariable() = default;
        
        virtual ~RandomVariable(){}
        
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            return this->operator()( dynamic_cast<const Sampler_T&>(S) );
        }
        
        virtual Real operator()( const Sampler_T & S ) const = 0;
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            return this->MinValue( dynamic_cast<const Sampler_T&>(S) );
        }
        
        virtual Real MinValue( const Sampler_T & S ) const = 0;
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            return this->MinValue( dynamic_cast<const Sampler_T&>(S) );
        }
        
        virtual Real MaxValue( const Sampler_T & S ) const = 0;

        
        __ADD_CLONE_CODE_FOR_ABSTRACT_CLASS__(RandomVariable)

        virtual RandomVariable & DownCast() override
        {
            return *this;
        }
        
        virtual const RandomVariable & DownCast() const override
        {
            return *this;
        }
        
    public:
        
        Int AmbientDimension() const override
        {
            return AmbDim;
        }
        
        virtual std::string Tag() const override = 0;
    };
    
} // namespace CycleSampler
