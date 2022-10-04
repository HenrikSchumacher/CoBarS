#pragma  once

namespace CyclicSampler {
    
#define CLASS RandomVariableBase

    template<
        typename Real    = double,
        typename Int     = long long
    >
    class CLASS
    {
        ASSERT_INT(Int);
        ASSERT_FLOAT(Real);
    
    public:
        
        using CyclicSampler_T = CyclicSamplerBase<Real,Int>;
        
        virtual ~CLASS() {}
            
        virtual Real operator()( const CyclicSampler_T & C ) const = 0;
        
        virtual Real MinValue( const CyclicSampler_T & C ) const = 0;
        
        virtual Real MaxValue( const CyclicSampler_T & C ) const = 0;
        
        CYCLICSAMPLER__ADD_CLONE_CODE_FOR_BASE_CLASS(CLASS)
        
    public:
        
        virtual CLASS & DownCast() = 0;

        virtual const CLASS & DownCast() const = 0;
        
        virtual Int AmbientDimension() const = 0;
        
        virtual bool RequiresSpaceCurve() const = 0;
        
        virtual std::string Tag() const = 0;
        
        virtual std::string ClassName() const 
        {
            return TO_STD_STRING(CLASS)+"<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
    };
#undef CLASS
        
} // namespace CyclicSampler
