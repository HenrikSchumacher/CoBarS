#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class GyradiusP;
    
    /*!
     * @brief Computes the p-gyradius of an instance of `CoBarS::SamplerBase<AmbDim,Real,Int>`. This is very expensive as it takes the l^p-norm of all chords of the polygon!
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class GyradiusP<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        GyradiusP( const Real exponent_ )
        :   exponent( exponent_ )
        {}
        
        // Copy constructor
        GyradiusP( const GyradiusP & other )
        :   exponent(other.exponent)
        {}
        
        // Move constructor
        GyradiusP( GyradiusP && other ) noexcept
        :
        exponent(other.exponent)
        {}
        
        virtual ~GyradiusP() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<GyradiusP> Clone () const
        {
            return std::shared_ptr<GyradiusP>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual GyradiusP * CloneImplementation() const override
        {
            return new GyradiusP(*this);
        }
        
    protected:
        
        const Real exponent = 2;
        
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            Real sum = 0;
            Real r2  = 0;
            
            const Real power = Frac<Real>(exponent,2);
            
            const Int n      = S.EdgeCount();
            
            for( Int k = 0; k < n; ++k )
            {
                Vector_T u = S.VertexPosition(k);
                
                for( Int l = k+1; l < n; ++l )
                {
                    Vector_T v = u;
                    
                    v -= S.VertexPosition(l);
                    
                    sum+= std::pow(v.SquaredNorm(),power);
                }
            }
            
            return std::pow( sum / (n * n), Inv<Real>(exponent) );
        }
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            (void)S;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            return Total(S.EdgeLengths()) * std::pow( S.EdgeCount(), -Inv<Real>(exponent) );
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("GyradiusP")+"("+ToString(exponent)+")";
        }
    };
    
} // namespace CoBarS
