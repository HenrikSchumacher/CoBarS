#pragma once

namespace CoBarS
{
        
#define CLASS BarycenterNorm
        
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
            // We treat the edges as massless.
            // All mass is concentrated in the vertices, and each vertex carries the same mass.
            
            const Weights_T & r = C.EdgeLengths();
            
            const Int n = C.EdgeCount();
            
            Vector_T b;
            
            b.SetZero();
            
            for( Int k = 0; k < n; ++ k )
            {
                for( Int i = 0; i < AmbDim; ++i )
                {
                    b[i] += r[i] * ( C.SpaceCoordinates(k,i) + C.SpaceCoordinates(k+1,i) );
                }
            }
            
            const Real factor = Scalar::Half<Real> / n;
            
            return b.Norm() * factor;
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return Scalar::Zero<Real>;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return Total(C.EdgeLengths());
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS);
        }
    };
        
#undef CLASS
        
} // namespace CoBarS
