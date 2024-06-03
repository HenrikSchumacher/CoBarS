#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class IterationCount;
    
    /*!
     * @brief Returns the number of iterations that `CoBarS::SamplerBase<AmbDim,Real,Int>` needed to converge.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class IterationCount<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        
        IterationCount() = default;
        
        virtual ~IterationCount() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<IterationCount> Clone () const
        {
            return std::shared_ptr<IterationCount>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual IterationCount * CloneImplementation() const override
        {
            return new IterationCount(*this);
        }
        
    protected:
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            return C.IterationCount();
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return C.MaxIterationCount();
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("IterationCount");
        }
    };
    
}  // namespace CoBarS
