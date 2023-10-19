#pragma once

namespace CoBarS
{
    
#define CLASS TotalCurvature
    
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
            
            const Int n              = C.EdgeCount();
            
            Real sum;
            
            // Handle wrap-around.
            {
                Vector_T u = C.EdgeCoordinates( n-1 );
                Vector_T v = C.EdgeCoordinates( 0   );
                
                const Real phi = AngleBetweenUnitVectors( u, v );
                
                sum = phi;
            }
            
            for( Int k = 0; k < n-1; ++k )
            {
                Vector_T u = C.EdgeCoordinates( k   );
                Vector_T v = C.EdgeCoordinates( k+1 );
                
                const Real phi = AngleBetweenUnitVectors( u, v );
                
                sum += phi;
            }
            
            return sum;
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return Scalar::Zero<Real>;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return C.EdgeCount() * Scalar::Pi<Real>;
        }
        
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
    
}  // namespace CoBarS
