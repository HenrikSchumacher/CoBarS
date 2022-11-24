#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class DiagonalLength : public RandomVariable<AmbDim,Real,Int>
{
public:
    
    using Sampler_T = Sampler<AmbDim,Real,Int>;
    
    DiagonalLength() = default;
    
    virtual ~DiagonalLength() override = default;
    
    __ADD_CLONE_CODE__(DiagonalLength)

protected:
    
    virtual Real operator()( const Sampler_T & C ) const override
    {
        Real r2 = static_cast<Real>(0);
        
        const Real * restrict const p = C.SpaceCoordinates();
        
        const Int last_vertex = C.EdgeCount()/2;
        
        for( Int i = 0; i < AmbDim; ++i )
        {
            const Real delta = p[AmbDim * last_vertex+i] - p[AmbDim*0+i];
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
        const Real * restrict const r = C.EdgeLengths();
        
        Real L = static_cast<Real>(0);
        
        const Int last_vertex = C.EdgeCount()/2;
        
        for( Int k = 0; k < last_vertex; ++k )
        {
            L += r[k];
        }
        
        return L;
    }
    
public:
    
    virtual std::string Tag() const  override
    {
        return "DiagonalLength";
    }
};
