#pragma once

namespace CoBarS
{
        
    template<typename SamplerBase_T> class BarycenterNorm;
    
    /*!
     * @brief Computes the Euclidean norm of the Euclidean barycenter of an instance of `CoBarS::SamplerBase<AmbDim,Real,Int>`.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class BarycenterNorm<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        BarycenterNorm() = default;
        
        virtual ~BarycenterNorm() override = default;

    public:
        
        [[nodiscard]] std::shared_ptr<BarycenterNorm> Clone () const
        {
            return std::shared_ptr<BarycenterNorm>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual BarycenterNorm * CloneImplementation() const override
        {
            return new BarycenterNorm(*this);
        }
        
    protected:
        
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            // We treat the edges as massless.
            // All mass is concentrated in the vertices, and each vertex carries the same mass.
            
            const Weights_T & r = S.EdgeLengths();
            
            const Int n = S.EdgeCount();
            
            Vector_T b;
            
            b.SetZero();
            
            for( Int i = 0; i < n; ++ i )
            {
                Vector_T u = S.VertexPosition(i  );
                Vector_T v = S.VertexPosition(i+1);
                
                const Real r_i = r[i];
                
                for( Int j = 0; j < AmbDim; ++j )
                {
                    b[j] += r_i * ( u[j] + v[j] );
                }
            }
            
            const Real factor = Scalar::Half<Real> / n;
            
            return b.Norm() * factor;
        }
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            (void)S;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            return Total(S.EdgeLengths());
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("BarycenterNorm");
        }
    };
        
} // namespace CoBarS
