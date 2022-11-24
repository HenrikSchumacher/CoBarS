#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class EdgeSpaceSamplingWeight : public RandomVariable<AmbDim,Real,Int>
{
public:
    
<<<<<<< HEAD
    using Sampler_T   = typename BASE::Sampler_T;
=======
    using Sampler_T = Sampler<AmbDim,Real,Int>;
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    
    EdgeSpaceSamplingWeight() = default;
    
    virtual ~EdgeSpaceSamplingWeight() override = default;
    
    __ADD_CLONE_CODE__(EdgeSpaceSamplingWeight)

protected:
    
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        return C.EdgeSpaceSamplingWeight();
    }
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
//            return static_cast<Real>(1)/( std::pow( C.EdgeCount(), AmbDim-1) );
        return static_cast<Real>(1)/( C.EdgeCount() );
    }
    
public:

    virtual std::string Tag() const  override
    {
        return "EdgeSpaceSamplingWeight";
    }
};
