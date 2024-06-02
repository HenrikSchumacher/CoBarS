#pragma once


namespace CoBarS
{
    template<typename SamplerBase_T> class EdgeQuotientSpaceSamplingWeight;
    
    template<int AmbDim, typename Real, typename Int>
    class EdgeQuotientSpaceSamplingWeight<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T        = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T     = typename Base_T::Weights_T;
        
        EdgeQuotientSpaceSamplingWeight() = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<EdgeQuotientSpaceSamplingWeight> Clone () const
        {
            return std::shared_ptr<EdgeQuotientSpaceSamplingWeight>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual EdgeQuotientSpaceSamplingWeight * CloneImplementation() const override
        {
            return new EdgeQuotientSpaceSamplingWeight(*this);
        }
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            return C.EdgeQuotientSpaceSamplingWeight();
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return Inv<Real>( C.EdgeCount() );
        }
        
    public:
        virtual std::string Tag() const  override
        {
            return std::string("EdgeQuotientSpaceSamplingWeight");
        }
    };
    
}  // namespace CoBarS
