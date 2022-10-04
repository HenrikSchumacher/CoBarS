#pragma once

#define CLASS TotalCurvature
#define BASE  RandomVariable<AmbDim,Real,Int>

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public BASE
{
public:
    
    using CyclicSampler_T   = typename BASE::CyclicSampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    using Weights_T         = typename BASE::Weights_T;
    
    CLASS() = default;
    
    virtual ~CLASS() override = default;
    
    __ADD_CLONE_CODE__(CLASS)

protected:
    
    virtual Real operator()( const CyclicSampler_T & C ) const override
    {
        
        const Int edge_count     = C.EdgeCount();
        const SpherePoints_T & y = C.EdgeCoordinates();

        Real sum;
        
        {
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( y.data(edge_count-1), y.data(0) );
            
            sum = phi;
        }
        
        for( Int k = 0; k < edge_count-1; ++k )
        {
            const Real phi = MyMath::AngleBetweenUnitVectors<AmbDim>( y.data(k), y.data(k+1) );
            
            sum += phi;
        }
        
        return sum;
    }
    
    virtual Real MinValue( const CyclicSampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const CyclicSampler_T & C ) const override
    {
        return C.EdgeCount() * static_cast<Real>(M_PI);
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return TO_STD_STRING(CLASS);
    }
    
    virtual std::string ClassName() const override
    {
        return TO_STD_STRING(CLASS)+"<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
    }
};

#undef BASE
#undef CLASS
