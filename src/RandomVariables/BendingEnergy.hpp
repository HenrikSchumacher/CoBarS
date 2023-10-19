#pragma once


namespace CoBarS
{
#define CLASS BendingEnergy
    
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
        
        explicit CLASS( const Real p_ )
        :   p( p_ )
        {}
        
        // Copy constructor
        explicit CLASS( const CLASS & other )
        :   p( other.p )
        {}
        
        // Move constructor
        explicit CLASS( CLASS && other ) noexcept
        :   p( other.p )
        {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        const Real p;
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            const Int n              = C.EdgeCount();
            const Weights_T & r      = C.EdgeLengths();
            
            Real sum;
            
            {
                const Real len = Scalar::Half<Real> * (r[n-1]+r[0]);
                
                Vector_T u = C.EdgeCoordinates(n-1);
                Vector_T v = C.EdgeCoordinates(0);
                
                const Real phi = AngleBetweenUnitVectors( u, v );
                
                sum = std::pow( phi / len, p ) * len;
            }
            
            for( Int k = 0; k < n-1; ++k )
            {
                const Real len = Scalar::Half<Real> * (r[k]+r[k+1]);
                
                Vector_T u = C.EdgeCoordinates(k);
                Vector_T v = C.EdgeCoordinates(k+1);
                
                const Real phi = AngleBetweenUnitVectors( u, v );
                
                sum += std::pow( phi / len, p ) * len;
            }
            
            return sum/p;
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return Scalar::Zero<Real>;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            const Int n         = C.EdgeCount();
            const Weights_T & r = C.EdgeLengths();
            
            Real sum;
            
            {
                const Real len = Scalar::Half<Real> * (r[n-1]+r[0]);
                
                const Real phi = Scalar::Pi<Real>;
                
                sum = std::pow( phi / len, p ) * len;
                
            }
            
            for( Int k = 0; k < n-1; ++k )
            {
                const Real len = Scalar::Half<Real> * (r[k]+r[k+1]);
                
                const Real phi = Scalar::Pi<Real>;
                
                sum += std::pow( phi / len, p ) * len;
            }
            
            return sum/p;
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS)+"("+ToString(p)+")";
        }
    };
    
#undef CLASS
    
} // namespace CoBarS
