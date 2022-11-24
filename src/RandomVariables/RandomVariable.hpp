#pragma once

<<<<<<< HEAD
namespace CycleSampler {

#define CLASS RandomVariable
#define BASE  RandomVariableBase<Real,Int>

=======
namespace CycleSampler
{
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    template<int AmbDim, typename Real = double, typename Int = long long>
    class RandomVariable : public RandomVariableBase<Real,Int>
    {
    public:
        
        using SamplerBase_T     = SamplerBase<Real,Int>;
        using Sampler_T         = Sampler<AmbDim,Real,Int>;
<<<<<<< HEAD
        using SpherePoints_T    = typename Sampler_T::SpherePoints_T;
        using SpacePoints_T     = typename Sampler_T::SpacePoints_T;
        using Weights_T         = typename Sampler_T::Weights_T;
        using Vector_T          = typename Sampler_T::Vector_T;
=======
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
        
        RandomVariable() = default;
        
        virtual ~RandomVariable(){}
        
<<<<<<< HEAD
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
=======
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            return this->operator()( dynamic_cast<const Sampler_T&>(S) );
        }
        
        virtual Real operator()( const Sampler_T & S ) const = 0;
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            return this->MinValue( dynamic_cast<const Sampler_T&>(S) );
        }
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
        
        virtual Real MinValue( const Sampler_T & S ) const = 0;
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            return this->MinValue( dynamic_cast<const Sampler_T&>(S) );
        }
        
        virtual Real MaxValue( const Sampler_T & S ) const = 0;

<<<<<<< HEAD
=======
        
        __ADD_CLONE_CODE_FOR_ABSTRACT_CLASS__(RandomVariable)
        
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    public:
        
        Int AmbientDimension() const override
        {
            return AmbDim;
        }
        
        virtual std::string Tag() const override = 0;
    };
    
} // namespace CycleSampler
