#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class ChordLength : public RandomVariable<AmbDim,Real,Int>
{
public:
    
<<<<<<< HEAD
    using Sampler_T         = typename BASE::Sampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    using Weights_T         = typename BASE::Weights_T;
=======
    using Sampler_T = Sampler<AmbDim,Real,Int>;
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    
    ChordLength( const Int first_vertex_, const Int last_vertex_)
    :   first_vertex( std::max( static_cast<Int>(0), first_vertex_) )
    ,   last_vertex(last_vertex_)
    {}
    
    // Copy constructor
    ChordLength( const ChordLength & other )
    :   first_vertex(other.first_vertex),
        last_vertex(other.last_vertex)
    {}

    // Move constructor
    ChordLength( ChordLength && other ) noexcept
    :
        first_vertex(other.first_vertex),
        last_vertex(other.last_vertex)
    {}
    
    virtual ~ChordLength() override = default;
    
    __ADD_CLONE_CODE__(ChordLength)

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
        
        const Real * restrict const p = C.SpaceCoordinates();
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            const Real delta = p[AmbDim*last_vertex+i] - p[AmbDim*first_vertex+i];
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
        const Real * restrict r = C.EdgeLengths();
        
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
        return "ChordLength("+ToString(first_vertex)+","+ToString(last_vertex)+")";
    }
};
