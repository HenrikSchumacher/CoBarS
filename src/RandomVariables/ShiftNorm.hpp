#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class ShiftNorm;
    
    /*!
     * @brief Computes the Euclidean norm of the conformal barycenter of an instance of `CoBarS::SamplerBase<AmbDim,Real,Int>`.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
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
        
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            return S.ShiftVector().Norm();
        }
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            (void)S;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            (void)S;
            
            return 1;
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("ShiftNorm");
        }
    };
    
}  // namespace CoBarS
