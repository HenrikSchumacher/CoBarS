#pragma once


namespace CoBarS
{
    template<typename SamplerBase_T> class DiagonalLength;
    
    /*!
     * @brief Computes a length of the chord between the first vertex and the middle vertex of an instance of `CoBarS::SamplerBase<AmbDim,Real,Int>`.
     *
     * @tparam AmbDim The dimension of the ambient space.
     *
     * @tparam Real A real floating point type.
     *
     * @tparam Int  An integer type.
     */
    
    template<int AmbDim, typename Real, typename Int>
    class DiagonalLength<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        DiagonalLength() = default;
        
        virtual ~DiagonalLength() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<DiagonalLength> Clone () const
        {
            return std::shared_ptr<DiagonalLength>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual DiagonalLength * CloneImplementation() const override
        {
            return new DiagonalLength(*this);
        }
        
    protected:
        
        virtual Real operator()( const SamplerBase_T & S ) const override
        {
            const Int last_vertex = S.EdgeCount()/2;
            
            Vector_T u = S.VertexPosition( last_vertex );
            Vector_T v = S.VertexPosition( 0           );
            
            u -= v;
            
            return u.Norm();
        }
        
        virtual Real MinValue( const SamplerBase_T & S ) const override
        {
            (void)S;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & S ) const override
        {
            const Weights_T & r = S.EdgeLengths();
            
            Real L = 0;
            
            const Int last_vertex = S.EdgeCount()/2;
            
            for( Int k = 0; k < last_vertex; ++k )
            {
                L += r[k];
            }
            
            return L;
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return "DiagonalLength";
        }
    };
    
}  // namespace CoBarS
