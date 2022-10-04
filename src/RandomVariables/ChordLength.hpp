#pragma  once

namespace CyclicSampler {

#define CLASS ChordLength
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
        
        CLASS( const Int first_vertex_, const Int last_vertex_)
        :   first_vertex( std::max( static_cast<Int>(0), first_vertex_) )
        ,   last_vertex(last_vertex_)
        {}
        
        // Copy constructor
        CLASS( const CLASS & other )
        :   first_vertex(other.first_vertex),
            last_vertex(other.last_vertex)
        {}

        // Move constructor
        CLASS( CLASS && other ) noexcept
        :
            first_vertex(other.first_vertex),
            last_vertex(other.last_vertex)
        {}
        
        virtual ~CLASS() override = default;
        
        __ADD_CLONE_CODE__(CLASS)

    protected:
        
        const Int first_vertex = 0;
        const Int last_vertex  = 0;
        
        virtual Real operator()( const CyclicSampler_T & C ) const override
        {
            Real r2 = static_cast<Real>(0);
            
            if( last_vertex > C.EdgeCount() )
            {
                return static_cast<Real>(0);
            }
            
            const SpacePoints_T & p = C.SpaceCoordinates();
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                const Real delta = p(last_vertex,i) - p(first_vertex,i);
                r2 += delta * delta;
            }
            
            return std::sqrt(r2);
        }
        
        virtual Real MinValue( const CyclicSampler_T & C ) const override
        {
            return static_cast<Real>(0);
        }
        
        virtual Real MaxValue( const CyclicSampler_T & C ) const override
        {
            auto & omega = C.Omega();
            
            Real L = static_cast<Real>(0);
            
            for( Int k = first_vertex; k < last_vertex; ++k )
            {
                L += omega[k];
            }
            
            return L;
        }
        
    public:
        
        virtual bool RequiresSpaceCurve() const override
        {
            return true;
        };
        
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS)+"("+ToString(first_vertex)+","+ToString(last_vertex)+")";
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
    };
        
    template<typename Real, typename Int>
    std::unique_ptr<ROOT> CONCAT(Make_,CLASS) (
        const Int amb_dim,
        const Int first_vertex_,
        const Int last_vertex_
    )
    {
        ROOT * r = nullptr;

        switch(amb_dim)
        {
            case 2:
            {
                r = new CLASS<2,Real,Int>(first_vertex_, last_vertex_);
                break;
            }
                
            case 3:
            {
                r = new CLASS<3,Real,Int>(first_vertex_, last_vertex_);
                break;
            }
                
            case 4:
            {
                r = new CLASS<4,Real,Int>(first_vertex_, last_vertex_);
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
