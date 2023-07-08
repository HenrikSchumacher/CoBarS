#pragma once

namespace CycleSampler
{
    
#define CLASS MaxAngle
    
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
            
            Real max_angle = static_cast<Real>(0);
            
            {
                Vector_T u = C.EdgeCoordinates( n-1 );
                Vector_T v = C.EdgeCoordinates( 0 );
                
                max_angle = std::max(max_angle, AngleBetweenUnitVectors( u, v ) );
            }
            
            for( Int k = 0; k < n-1; ++k )
            {
                Vector_T u = C.EdgeCoordinates( k );
                Vector_T v = C.EdgeCoordinates( k+1 );
                
                max_angle = std::max(max_angle, AngleBetweenUnitVectors( u, v ) );
            }
            
            return max_angle;
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return static_cast<Real>(0);
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return Scalar::Pi<Real>;
        }
        
    public:
        
        virtual std::string Tag() const override
        {
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
}
