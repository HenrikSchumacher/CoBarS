#pragma once


namespace CycleSampler
{
    
#define CLASS ChordLength
    
    template<typename SamplerBase_T> class CLASS;
    
    template<int AmbDim, typename Real, typename Int>
    class CLASS<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
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
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
             if( last_vertex > C.EdgeCount() )
            {
                return Scalar::Zero<Real>;
            }
            
            Vector_T u = C.SpaceCoordinates( last_vertex  );
            Vector_T v = C.SpaceCoordinates( first_vertex );
            
            u -= v;
            
            return u.Norm();
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            return Scalar::Zero<Real>;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            const Weights_T & r = C.EdgeLengths();
            
            Real L = Scalar::Zero<Real>;
            
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
    
#undef CLASS
    
}
