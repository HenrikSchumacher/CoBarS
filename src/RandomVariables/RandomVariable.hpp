#pragma  once

namespace CyclicSampler {

#define CLASS RandomVariable
#define BASE  RandomVariableBase<Real,Int>

    template<int AmbDim, typename Real = double, typename Int = long long>
    class CLASS : public BASE
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using CyclicSampler_T   = typename BASE::CyclicSampler_T;
        using SpherePoints_T    = typename CyclicSampler_T::SpherePoints_T;
        using SpacePoints_T     = typename CyclicSampler_T::SpacePoints_T;
        using Weights_T         = typename CyclicSampler_T::Weights_T;
        
        CLASS() = default;
        
        virtual ~CLASS(){}
        
        virtual Real operator()( const CyclicSampler_T & C ) const override = 0;
        
        virtual Real MinValue( const CyclicSampler_T & C ) const override = 0;
        
        virtual Real MaxValue( const CyclicSampler_T & C ) const override = 0;
        
        __ADD_CLONE_CODE_FOR_ABSTRACT_CLASS__(CLASS)

        virtual CLASS & DownCast() override
        {
            return *this;
        }
        
        virtual const CLASS & DownCast() const override
        {
            return *this;
        }
        
    public:
        
        Int AmbientDimension() const override
        {
            return AmbDim;
        }
        
        virtual bool RequiresSpaceCurve() const override = 0;
        
        virtual std::string Tag() const override = 0;
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
    };
        
#undef BASE
#undef CLASS
    
} // namespace CyclicSampler
