#pragma once

namespace CycleSampler {
    
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
        
        using SamplerBase_T = SamplerBase<Real,Int>;
        
        virtual ~CLASS() {}
            
        virtual Real operator()( const SamplerBase_T & C ) const = 0;
        
        virtual Real MinValue( const SamplerBase_T & C ) const = 0;
        
        virtual Real MaxValue( const SamplerBase_T & C ) const = 0;
        
        __ADD_CLONE_CODE_FOR_BASE_CLASS__(CLASS)
        
    public:
        
        virtual Int AmbientDimension() const = 0;
        
        virtual std::string Tag() const = 0;
    };
#undef CLASS
        
} // namespace CycleSampler
