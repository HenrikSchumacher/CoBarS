#pragma once

namespace CycleSampler
{
    template< typename Real = double, typename Int = long long >
    class RandomVariableBase
    {
        ASSERT_INT(Int);
        ASSERT_FLOAT(Real);
    
    public:
        
        using Sampler_T = SamplerBase<Real,Int>;
        
        virtual ~RandomVariableBase() {}
            
        virtual Real operator()( const Sampler_T & C ) const = 0;
        
        virtual Real MinValue( const Sampler_T & C ) const = 0;
        
        virtual Real MaxValue( const Sampler_T & C ) const = 0;
        
        __ADD_CLONE_CODE_FOR_BASE_CLASS__(RandomVariableBase)
        
    public:
        
        virtual Int AmbientDimension() const = 0;
        
        virtual std::string Tag() const = 0;
    };
        
} // namespace CycleSampler
