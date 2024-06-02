#pragma once

namespace CoBarS
{
    template<typename SamplerBase_T> class Gyradius;
    
    template<int AmbDim, typename Real, typename Int>
    class Gyradius<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        
        Gyradius() = default;
        
        virtual ~Gyradius() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<Gyradius> Clone () const
        {
            return std::shared_ptr<Gyradius>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual Gyradius * CloneImplementation() const override
        {
            return new Gyradius(*this);
        }
        
    protected:
        
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            Real r2 = 0;
            
            const Int n = C.EdgeCount();

            for( Int k = 0; k < n; ++k )
            {
                r2 += C.SpaceCoordinates(k).SquaredNorm();
            }
            
            return std::sqrt( r2/n );
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            return Total( C.EdgeLengths() ) / std::sqrt( static_cast<Real>(C.EdgeCount()) );
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("Gyradius");
        }
    };
    
} // namespace CoBarS
