#pragma once


namespace CoBarS
{
    template<typename SamplerBase_T> class EdgeSpaceSamplingWeight;
    
    template<int AmbDim, typename Real, typename Int>
    class EdgeSpaceSamplingWeight<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        
        EdgeSpaceSamplingWeight() = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<EdgeSpaceSamplingWeight> Clone () const
        {
            return std::shared_ptr<EdgeSpaceSamplingWeight>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual EdgeSpaceSamplingWeight * CloneImplementation() const override
        {
            return new EdgeSpaceSamplingWeight(*this);
        }
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            return C.EdgeSpaceSamplingWeight();
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
            return std::string("EdgeSpaceSamplingWeight");
        }
    };
    
}  // namespace CoBarS
