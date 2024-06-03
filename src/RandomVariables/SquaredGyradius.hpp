#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class SquaredGyradius;
    
    /*!
     * @brief Computes the squared gyradius of an instance of `CoBarS::SamplerBase<AmbDim,Real,Int>`.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class SquaredGyradius<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        
        SquaredGyradius() = default;
        
        virtual ~SquaredGyradius() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<SquaredGyradius> Clone () const
        {
            return std::shared_ptr<SquaredGyradius>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual SquaredGyradius * CloneImplementation() const override
        {
            return new SquaredGyradius(*this);
        }
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            Real r2 = 0;
            
            const Int n = C.EdgeCount();

            for( Int k = 0; k < n; ++k )
            {
                r2 += C.VertexPosition(k).SquaredNorm();
            }
            
            return r2/n;
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            const Real L = Total( C.EdgeLengths() );
            
            return L * L / static_cast<Real>(C.EdgeCount());
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("SquaredGyradius");
        }
    };
}  // namespace CoBarS

