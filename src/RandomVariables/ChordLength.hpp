#pragma once

#define CLASS ChordLength

template<int AmbDim, typename Real = double, typename Int = long long>
class CLASS : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Base_T            = RandomVariable<AmbDim,Real,Int>;
    using Sampler_T         = typename Base_T::Sampler_T;
    using SpherePoints_T    = typename Base_T::SpherePoints_T;
    using SpacePoints_T     = typename Base_T::SpacePoints_T;
    using Weights_T         = typename Base_T::Weights_T;
    
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
    
    virtual Real operator()( const Sampler_T & C ) const override
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
    
    virtual Real MinValue( const Sampler_T & C ) const override
    {
        return static_cast<Real>(0);
    }
    
    virtual Real MaxValue( const Sampler_T & C ) const override
    {
        const Weights_T & r = C.EdgeLengths();
        
        Real L = static_cast<Real>(0);
        
        for( Int k = first_vertex; k < last_vertex; ++k )
        {
            L += r[k];
        }
        
        return L;
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return TO_STD_STRING(CLASS)+"("+ToString(first_vertex)+","+ToString(last_vertex)+")";
    }
};
        
#undef BASE
#undef CLASS
