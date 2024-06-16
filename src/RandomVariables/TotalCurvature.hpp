#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class TotalCurvature;
    
    /*!
     * @brief Computes the total curvature of an instance of `CoBarS::SamplerBase<AmbDim,Real,Int>`.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class TotalCurvature<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        TotalCurvature() = default;
        
        virtual ~TotalCurvature() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<TotalCurvature> Clone () const
        {
            return std::shared_ptr<TotalCurvature>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual TotalCurvature * CloneImplementation() const override
        {
            return new TotalCurvature(*this);
        }
        
    protected:
        
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            
            const Int n = S.EdgeCount();
            
            Real sum;
            
            // Handle wrap-around.
            {
                const Real phi = AngleBetweenUnitVectors(
                    S.EdgeVector(n-1), S.EdgeVector(0)
                );
                
                sum = phi;
            }
            
            for( Int k = 0; k < n-1; ++k )
            {
                const Real phi = AngleBetweenUnitVectors(
                    S.EdgeVector(k), S.EdgeVector(k+1)
                );
                
                sum += phi;
            }
            
            return sum;
        }
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            (void)S;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            return S.EdgeCount() * Scalar::Pi<Real>;
        }
        
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("TotalCurvature");
        }
    };
    
    
}  // namespace CoBarS
