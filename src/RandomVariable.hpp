#pragma once

namespace CoBarS
{
    /*!
     * @brief The base class for all random variables of `CoBarS::SamplerBase<AmbDim,Real,Int>`. This enables runtime polymorphism, so that template parameters  of `CoBarS::Sampler` like the pseudorandom number generator do not matter.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
        static_assert(FloatQ<Real>,"");
        static_assert(IntQ<Int>,"");
        
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        using Weights_T         = typename SamplerBase_T::Weights_T;
        using Vector_T          = typename SamplerBase_T::Vector_T;
        
        RandomVariable() = default;
        
        virtual ~RandomVariable(){}
        
        virtual Real operator()( const SamplerBase_T & C ) const = 0;
        
        virtual Real MinValue( const SamplerBase_T & C ) const = 0;
        
        virtual Real MaxValue( const SamplerBase_T & C ) const = 0;
        
    public:
        
        [[nodiscard]] std::shared_ptr<RandomVariable> Clone () const
        {
            return std::shared_ptr<RandomVariable>(CloneImplementation());
        }
        
    private:
        
        [[nodiscard]] virtual RandomVariable * CloneImplementation() const = 0;

    public:
        
        Int AmbientDimension() const
        {
            return AmbDim;
        }
        
        virtual std::string Tag() const = 0;
        
    }; // RandomVariable
    
} // namespace CoBarS
