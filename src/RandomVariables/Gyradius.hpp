#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class Gyradius;
    
    /*!
     * @brief Computes the gyradius of an instance of `CoBarS::SamplerBase<AmbDim,Real,Int>`.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class Gyradius<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        
        Gyradius() = default;
        
        virtual ~Gyradius() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<Gyradius> Clone() const
        {
            return std::shared_ptr<Gyradius>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual Gyradius * CloneImplementation() const override
        {
            return new Gyradius(*this);
        }
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            Real r2 = 0;
            
            const Int n = S.EdgeCount();

            for( Int k = 0; k < n; ++k )
            {
                r2 += S.VertexPosition(k).SquaredNorm();
            }
            
            return std::sqrt( r2/n );
        }
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            (void)S;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            return Total( S.EdgeLengths() ) / std::sqrt( static_cast<Real>(S.EdgeCount()) );
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("Gyradius");
        }
    };
    
} // namespace CoBarS
