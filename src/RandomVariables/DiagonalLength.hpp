#pragma once


namespace CoBarS
{
    
#define CLASS DiagonalLength
    
    template<typename SamplerBase_T> class CLASS;
    
    template<int AmbDim, typename Real, typename Int>
    class CLASS<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        CLASS() = default;
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
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
            return Scalar::Zero<Real>;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            const Weights_T & r = C.EdgeLengths();
            
            Real L = Scalar::Zero<Real>;
            
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
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
    
}  // namespace CoBarS
