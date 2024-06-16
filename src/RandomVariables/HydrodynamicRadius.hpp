#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class HydrodynamicRadius;
    
    /*!
     * @brief Returns the hygrodynamic radius of of an instance of `CoBarS::SamplerBase<AmbDim,Real,Int>`.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class HydrodynamicRadius<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        HydrodynamicRadius() = default;
        
        virtual ~HydrodynamicRadius() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<HydrodynamicRadius> Clone () const
        {
            return std::shared_ptr<HydrodynamicRadius>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual HydrodynamicRadius * CloneImplementation() const override
        {
            return new HydrodynamicRadius(*this);
        }
        
        static constexpr Real eps = std::numeric_limits<Real>::min();
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            Real sum = 0;
            
            const Int n = S.EdgeCount();
            
            for( Int k = 0; k < n; ++k )
            {
                Vector_T u = S.VertexPosition(k);
                
                for( Int l = k+1; l < n; ++l )
                {
                    Vector_T v = u;
                    
                    v -= S.VertexPosition(l);
                                        
                    sum+= Inv<Real>(v.Norm() + eps);
                }
            }
            
            return (n * n)/sum;
        }
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            (void)S;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            return S.EdgeLengths().Total();
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("HydrodynamicRadius");
        }
    };
    
} // namespace CoBarS
