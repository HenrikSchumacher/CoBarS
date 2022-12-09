#pragma once

namespace CycleSampler
{
    
#define CLASS IterationCount
    
    template<int AmbDim, typename Real = double, typename Int = long long>
    class CLASS : public RandomVariable<AmbDim,Real,Int>
    {
    private:
        
        using Base_T    = RandomVariable<AmbDim,Real,Int>;
        
    public:
        
        using Sampler_T = typename Base_T::Sampler_T;
        
        CLASS() = default;
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        virtual Real operator()( const Sampler_T & C ) const override
        {
            return C.IterationCount();
        }
        
        virtual Real MinValue( const Sampler_T & C ) const override
        {
            return static_cast<Real>(0);
        }
        
        virtual Real MaxValue( const Sampler_T & C ) const override
        {
            return C.MaxIterationCount();
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
}
