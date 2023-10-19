#pragma once

namespace CoBarS
{
    template<int AmbDim, typename Real, typename Int>
    class RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        using Weights_T         = typename SamplerBase_T::Weights_T;
        using Vector_T          = typename SamplerBase_T::Vector_T;
        
        RandomVariable() = default;
        
        virtual ~RandomVariable(){}
        
        virtual Real operator()( const SamplerBase_T & C ) const = 0;
        
        virtual Real MinValue( const SamplerBase_T & C ) const = 0;
        
        virtual Real MaxValue( const SamplerBase_T & C ) const = 0;
        
        __ADD_CLONE_CODE_FOR_BASE_CLASS__(RandomVariable)

    public:
        
        Int AmbientDimension() const
        {
            return AmbDim;
        }
        
        virtual std::string Tag() const = 0;
        
    }; // RandomVariable
    
} // namespace CoBarS
