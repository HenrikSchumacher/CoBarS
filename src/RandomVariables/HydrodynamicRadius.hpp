#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class HydrodynamicRadius;
    
    template<int AmbDim, typename Real, typename Int>
    class HydrodynamicRadius<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        HydrodynamicRadius() = default;
        
        virtual ~HydrodynamicRadius() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<HydrodynamicRadius> Clone () const
        {
            return std::shared_ptr<HydrodynamicRadius>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual HydrodynamicRadius * CloneImplementation() const override
        {
            return new HydrodynamicRadius(*this);
        }
        
        static constexpr Real eps = std::numeric_limits<Real>::min();
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            Real sum = 0;
            Real r2  = 0;
            
            const Int n = C.EdgeCount();
            
            for( Int k = 0; k < n; ++k )
            {
                Vector_T u = C.SpaceCoordinates(k);
                
                for( Int l = k+1; l < n; ++l )
                {
                    r2 = 0;
                    
                    Vector_T v = u;
                    
                    v -= C.SpaceCoordinates(l);
                                        
                    sum+= Inv<Real>(v.Norm() + eps);
                }
            }
            
            return (n * n)/sum;
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return C.EdgeLengths().Total();
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("HydrodynamicRadius");
        }
    };
    
} // namespace CoBarS
