#pragma  once

namespace CyclicSampler {

#define CLASS BendingEnergy
#define BASE  RandomVariable<AmbDim,Real,Int>
#define ROOT  RandomVariableBase<Real,Int>

    template<int AmbDim, typename Real = double, typename Int = long long>
    class CLASS : public BASE
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using CyclicSampler_T   = typename BASE::CyclicSampler_T;
        using SpherePoints_T    = typename BASE::SpherePoints_T;
        using SpacePoints_T     = typename BASE::SpacePoints_T;
        using Weights_T         = typename BASE::Weights_T;
        
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
        
        CYCLICSAMPLER__ADD_CLONE_CODE(CLASS)

    protected:
        
        const Real p;
        
        virtual Real operator()( const CyclicSampler_T & C ) const override
        {
            const Int edge_count     = C.EdgeCount();
            const SpherePoints_T & y = C.EdgeCoordinates();
            const Weights_T & omega = C.Omega();
        
            Real sum;
            
            {
                const Real len = static_cast<Real>(0.5)*(omega[edge_count-1]+omega[0]);
                
                const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( y.data(edge_count-1), y.data(0) );
                
                sum = std::pow( phi / len, p ) * len;
            }
            
            for( Int k = 0; k < edge_count-1; ++k )
            {
                const Real len = static_cast<Real>(0.5)*(omega[k]+omega[k+1]);
                
                const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( y.data(k), y.data(k+1) );
                
                sum += std::pow( phi / len, p ) * len;
            }
            
            return sum/p;
        }
        
        virtual Real MinValue( const CyclicSampler_T & C ) const override
        {
            return static_cast<Real>(0);
        }
        
        virtual Real MaxValue( const CyclicSampler_T & C ) const override
        {
            const Int edge_count     = C.EdgeCount();
            const Weights_T & omega = C.Omega();
        
            Real sum;
            
            {
                const Real len = static_cast<Real>(0.5)*(omega[edge_count-1]+omega[0]);
                
                const Real phi = static_cast<Real>(M_PI);
                
                sum = std::pow( phi / len, p ) * len;
                
            }
            
            for( Int k = 0; k < edge_count-1; ++k )
            {
                const Real len = static_cast<Real>(0.5)*(omega[k]+omega[k+1]);
                
                const Real phi = static_cast<Real>(M_PI);
                
                sum += std::pow( phi / len, p ) * len;
            }
            
            return sum/p;
        }
        
    public:
        
        virtual bool RequiresSpaceCurve() const override
        {
            return false;
        };
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS)+"("+ToString(p)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
    };
        
    template<typename Real, typename Int>
    std::unique_ptr<ROOT> CONCAT(Make_,CLASS) (
        const Int amb_dim,
        const Real p = static_cast<Real>(2)
    )
    {
        ROOT * r = nullptr;

        switch(amb_dim)
        {
            case 2:
            {
                r = new CLASS<2,Real,Int>(p);
                break;
            }
                
            case 3:
            {
                r = new CLASS<3,Real,Int>(p);
                break;
            }
                
            case 4:
            {
                r = new CLASS<4,Real,Int>(p);
                break;
            }
            default:
            {
                eprint( "Make_"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported. Returning nullptr.");
                r = nullptr;
                break;
            }
        }
        
        return std::unique_ptr<ROOT>(r);
    }
    
#undef ROOT
#undef BASE
#undef CLASS
    
} // namespace CyclicSampler
