#pragma once

namespace CoBarS
{

#define CLASS Gyradius
    
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
        
        CLASS() = default;
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            Real r2 = Scalar::Zero<Real>;
            
            const Int n = C.EdgeCount();

            for( Int k = 0; k < n; ++k )
            {
                r2 += C.SpaceCoordinates(k).SquaredNorm();
            }
            
            return std::sqrt( r2/n );
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return Scalar::Zero<Real>;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return Total( C.EdgeLengths() ) / std::sqrt( static_cast<Real>(C.EdgeCount()) );
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
    
} // namespace CoBarS
