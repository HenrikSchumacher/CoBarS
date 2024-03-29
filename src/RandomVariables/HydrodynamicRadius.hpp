#pragma once

namespace CoBarS
{
    
#define CLASS HydrodynamicRadius
    
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
        
        static constexpr Real eps = std::numeric_limits<Real>::min();
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            Real sum = Scalar::Zero<Real>;
            Real r2  = Scalar::Zero<Real>;
            
            const Int n = C.EdgeCount();
            
            for( Int k = 0; k < n; ++k )
            {
                Vector_T u = C.SpaceCoordinates(k);
                
                for( Int l = k+1; l < n; ++l )
                {
                    r2 = Scalar::Zero<Real>;
                    
                    Vector_T v = u;
                    
                    v -= C.SpaceCoordinates(l);
                                        
                    sum+= Inv<Real>(v.Norm() + eps);
                }
            }
            
            return (n * n)/sum;
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return Scalar::Zero<Real>;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return C.EdgeLengths().Total();
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
    
} // namespace CoBarS
