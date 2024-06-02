#pragma once


namespace CoBarS
{
    template<typename SamplerBase_T> class EdgeQuotientSpaceSamplingCorrection;
    
    template<int AmbDim, typename Real, typename Int>
    class EdgeQuotientSpaceSamplingCorrection<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T        = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T     = typename Base_T::Weights_T;
        
        EdgeQuotientSpaceSamplingCorrection() = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<EdgeQuotientSpaceSamplingCorrection> Clone () const
        {
            return std::shared_ptr<EdgeQuotientSpaceSamplingCorrection>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual EdgeQuotientSpaceSamplingCorrection * CloneImplementation() const override
        {
            return new EdgeQuotientSpaceSamplingCorrection(*this);
        }
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            return C.EdgeQuotientSpaceSamplingCorrection();
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 2;
        }
        
    public:
        virtual std::string Tag() const  override
        {
            return std::string("EdgeQuotientSpaceSamplingCorrection");
        }
    };
    
}  // namespace CoBarS

