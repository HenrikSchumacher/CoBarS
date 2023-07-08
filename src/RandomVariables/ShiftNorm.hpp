#pragma once

namespace CycleSampler
{
    
#define CLASS ShiftNorm
    
    template<typename SamplerBase_T> class CLASS;
    
    template<int AmbDim, typename Real, typename Int>
    class CLASS<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        
        CLASS() = default;
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            return C.ShiftVector().Norm();
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return static_cast<Real>(0);
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return static_cast<Real>(1);
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
}
