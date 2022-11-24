#pragma once

namespace CycleSampler
{
    
#define CLASS RandomVariableBase

    template<
        typename Real    = double,
        typename Int     = long long
    >
    class CLASS
    {
        ASSERT_INT(Int);
        ASSERT_FLOAT(Real);
    
    public:
        
        using Sampler_T = SamplerBase<Real,Int>;
        
        virtual ~CLASS() {}
            
        virtual Real operator()( const Sampler_T & C ) const = 0;
        
        virtual Real MinValue( const Sampler_T & C ) const = 0;
        
        virtual Real MaxValue( const Sampler_T & C ) const = 0;
        
        __ADD_CLONE_CODE_FOR_BASE_CLASS__(CLASS)
        
    public:
        
        virtual CLASS & DownCast() = 0;

        virtual const CLASS & DownCast() const = 0;
        
        virtual Int AmbientDimension() const = 0;
        
        virtual std::string Tag() const = 0;
    };
#undef CLASS
        
} // namespace CycleSampler
