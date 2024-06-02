#pragma once


namespace CoBarS
{
    template<typename SamplerBase_T> class ChordLength;
    
    template<int AmbDim, typename Real, typename Int>
    class ChordLength<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        ChordLength( const Int first_vertex_, const Int last_vertex_)
        :   first_vertex( std::max( static_cast<Int>(0), first_vertex_) )
        ,   last_vertex ( std::max( static_cast<Int>(0), last_vertex_ ) )
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
        last_vertex (other.last_vertex )
        {}
        
        virtual ~ChordLength() override = default;
        
    public:
        
        [[nodiscard]] std::shared_ptr<ChordLength> Clone () const
        {
            return std::shared_ptr<ChordLength>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual ChordLength * CloneImplementation() const override
        {
            return new ChordLength(*this);
        }
        
    protected:
        
        const Int first_vertex = 0;
        const Int last_vertex  = 0;
        
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            if( last_vertex > C.EdgeCount() )
            {
                return 0;
            }
            
            Vector_T u = C.SpaceCoordinates( last_vertex  );
            Vector_T v = C.SpaceCoordinates( first_vertex );
            
            u -= v;
            
            return u.Norm();
        }
        
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 0;
        }
        
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            const Weights_T & r = C.EdgeLengths();
            
            Real L = 0;
            
            for( Int k = first_vertex; k < last_vertex; ++k )
            {
                L += r[k];
            }
            
            return L;
        }
        
    public:
        
        virtual std::string Tag() const  override
        {
            return std::string("ChordLength")+"("+ToString(first_vertex)+","+ToString(last_vertex)+")";
        }
    };
    
}  // namespace CoBarS
