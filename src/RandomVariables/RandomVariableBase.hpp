#pragma once

<<<<<<< HEAD
namespace CycleSampler {
    
#define CLASS RandomVariableBase

    template<
        typename Real    = double,
        typename Int     = long long
    >
    class CLASS
=======
namespace CycleSampler
{
    template< typename Real = double, typename Int = long long >
    class RandomVariableBase
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    {
        ASSERT_INT(Int);
        ASSERT_FLOAT(Real);
    
    public:
        
<<<<<<< HEAD
        using SamplerBase_T = SamplerBase<Real,Int>;
=======
        using Sampler_T = SamplerBase<Real,Int>;
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
        
        virtual ~RandomVariableBase() {}
            
<<<<<<< HEAD
        virtual Real operator()( const SamplerBase_T & C ) const = 0;
        
        virtual Real MinValue( const SamplerBase_T & C ) const = 0;
        
        virtual Real MaxValue( const SamplerBase_T & C ) const = 0;
=======
        virtual Real operator()( const Sampler_T & C ) const = 0;
        
        virtual Real MinValue( const Sampler_T & C ) const = 0;
        
        virtual Real MaxValue( const Sampler_T & C ) const = 0;
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
        
        __ADD_CLONE_CODE_FOR_BASE_CLASS__(RandomVariableBase)
        
    public:
        
        virtual Int AmbientDimension() const = 0;
        
        virtual std::string Tag() const = 0;
    };
        
} // namespace CycleSampler
