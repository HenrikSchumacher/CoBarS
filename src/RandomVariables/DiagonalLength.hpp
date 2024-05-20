#pragma once


namespace CoBarS
{
    template<typename SamplerBase_T> class DiagonalLength;
    
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
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            const Int last_vertex = C.EdgeCount()/2;
            
            Vector_T u = C.SpaceCoordinates( last_vertex );
            Vector_T v = C.SpaceCoordinates( 0           );
            
            u -= v;
            
            return u.Norm();
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            const Weights_T & r = C.EdgeLengths();
            
            Real L = 0;
            
            const Int last_vertex = C.EdgeCount()/2;
            
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
