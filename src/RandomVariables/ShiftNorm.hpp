#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class ShiftNorm;
    
    template<int AmbDim, typename Real, typename Int>
    class ShiftNorm<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        
        ShiftNorm() = default;
        
        virtual ~ShiftNorm() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<ShiftNorm> Clone () const
        {
            return std::shared_ptr<ShiftNorm>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual ShiftNorm * CloneImplementation() const override
        {
            return new ShiftNorm(*this);
        }
        
    protected:
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            return C.ShiftVector().Norm();
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return 1;
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("ShiftNorm");
        }
    };
    
}  // namespace CoBarS
