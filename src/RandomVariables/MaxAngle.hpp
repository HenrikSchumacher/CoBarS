#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class MaxAngle;
    
    template<int AmbDim, typename Real, typename Int>
    class MaxAngle<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        MaxAngle() = default;
        
        virtual ~MaxAngle() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<MaxAngle> Clone () const
        {
            return std::shared_ptr<MaxAngle>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual MaxAngle * CloneImplementation() const override
        {
            return new MaxAngle(*this);
        }
        
    protected:
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            
            const Int n    = C.EdgeCount();
            
            Real max_angle = 0;
            
            {
                Vector_T u = C.EdgeCoordinates( n-1 );
                Vector_T v = C.EdgeCoordinates( 0   );
                
                max_angle = std::max(max_angle, AngleBetweenUnitVectors( u, v ) );
            }
            
            for( Int k = 0; k < n-1; ++k )
            {
                Vector_T u = C.EdgeCoordinates( k   );
                Vector_T v = C.EdgeCoordinates( k+1 );
                
                max_angle = std::max(max_angle, AngleBetweenUnitVectors( u, v ) );
            }
            
            return max_angle;
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return Scalar::Pi<Real>;
        }
        
    public:
        
        virtual std::string Tag() const override
        {
            return std::string("MaxAngle");
        }
    };
    
}  // namespace CoBarS
