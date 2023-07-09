#pragma once

namespace CycleSampler
{
    
#define CLASS GyradiusP
    
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
        
        CLASS( const Real exponent_ )
        :   exponent( exponent_ )
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   exponent(other.exponent)
        {}
        
        // Move constructor
        CLASS( CLASS && other ) noexcept
        :
        exponent(other.exponent)
        {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        const Real exponent = 2;
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            Real sum = Scalar::Zero<Real>;
            Real r2  = Scalar::Zero<Real>;
            
            const Real power = Scalar::Half<Real> * exponent;
            
            const Int n      = C.EdgeCount();
            
            for( Int k = 0; k < n; ++k )
            {
                for( Int l = k+1; l < n; ++l )
                {
                    r2 = Scalar::Zero<Real>;
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        const Real delta = C.SpaceCoordinates(l,i) - C.SpaceCoordinates(k,i);
                        
                        r2 += delta * delta;
                    }
                    
                    sum+= std::pow(r2,power);
                }
            }
            
            return std::pow( sum / (n * n), Scalar::Inv<Real>(exponent) );
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return Scalar::Zero<Real>;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return Total(C.EdgeLengths()) * std::pow( C.EdgeCount(), -Scalar::Inv<Real>(exponent) );
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS)+"("+ToString(exponent)+")";
        }
    };
    
#undef CLASS
    
}
