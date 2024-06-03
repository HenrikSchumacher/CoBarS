#pragma once

namespace CoBarS
{
    
    template<
        int AmbDim_,
        typename Real_    = double,
        typename Int_     = int_fast32_t
    >
    class DouadyEarleExtension
    {
        // Load a polyline, project the vertices onto the sphere and allow for evaluation of the Douady-Earle extension
        
        static_assert(FloatQ<Real_>,"");
        static_assert(Scalar::RealQ<Real_>,"");
        static_assert(IntQ<Int_>,"");
        
    public:

        using Real   = Real_;
        using Int    = Int_;
        
        static constexpr Int AmbDim = int_cast<Int>(AmbDim_);
        
        // It does not matter which pseudorandom number generator we use.
        // Nothing is random here.
        using Prng_T  = Xoshiro256Plus;
        
        using Sampler2D_T        = Sampler<2     ,Real,Int,Prng_T,true,true>;
        using Sampler_T          = Sampler<AmbDim,Real,Int,Prng_T,true,true>;

        using Vector2D_T         = typename Sampler2D_T::Vector_T;
        using Vector_T           = typename   Sampler_T::Vector_T;
        
        using VectorList2D_T     = typename Sampler2D_T::VectorList_T;
        using VectorList_T       = typename   Sampler_T::VectorList_T;
        
        

//        using SquareMatrix_T     = typename Base_T::SquareMatrix_T;
//        using SymmetricMatrix_T  = typename Base_T::SymmetricMatrix_T;
//
//        using RandomVariable_T   = typename Base_T::RandomVariable_T;
//        using RandomVariable_Ptr = std::shared_ptr<RandomVariable_T>;
//
//        using Weights_T          = typename Base_T::Weights_T;
//        using Setting_T          = typename Base_T::Setting_T;
//

//        using Matrix_T           = Tensor2<Real,Int>;
    
    private:
        

        const Int edge_count_;
        
        mutable Sampler2D_T S2D_;
        mutable Sampler_T   S_;
        
        mutable Int vertex_count_;
        mutable VectorList_T curve_;
        
    public:
        
        DouadyEarleExtension() = default;
        
        ~DouadyEarleExtension() = default;
        
        explicit DouadyEarleExtension(
            const Int edge_count
        )
        :   edge_count_ ( edge_count  )
        ,   S2D_        ( edge_count_ )
        ,   S_          ( edge_count_ )
        {
            const Real step  = Scalar::TwoPi<Real> / edge_count_;
            
            for( Int k = 0; k < edge_count_; ++k )
            {
                const Real angle = step * k;
                
                Vector2D_T z = {std::cos(angle), std::sin(angle) };
                
                S2D_.ReadInitialEdgeVectors(z,k,false);
            }
        }
        
        
        // Copy constructor
        DouadyEarleExtension( const DouadyEarleExtension & other )
        :   edge_count_   ( other.edge_count_   )
        ,   S2D_          ( other.S2D_          )
        ,   S_            ( other.S_            )
        ,   vertex_count_ ( other.vertex_count_ )
        ,   curve_        ( other.curve_        )
        {}
        
        void LoadCurve( 
            cptr<Real> curve,
            const Int vertex_count,
            const bool normalizeQ = false
        )
        {
            vertex_count_ = vertex_count;
            
            if( curve_.Dimension(1) != vertex_count_ )
            {
                curve_ = VectorList_T( vertex_count_ );
            }
            
            if( normalizeQ )
            {
                for( Int k = 0; k < vertex_count_; ++k )
                {
                    Vector_T u;
                    u.Read( &curve[AmbDim * k] );
                    u.Normalize();
                    u.Write(curve_,k);
                }
            }
            else
            {
                for( Int k = 0; k < vertex_count_; ++k )
                {
                    Vector_T u;
                    u.Read( &curve[AmbDim * k] );
                    u.Write(curve_,k);
                }
            }
        }
        
        void operator()( cptr<Real> w_in, mptr<Real> w_out )
        {
            // Take 2D input vector and shift it to 0 so that
            Vector2D_T w;
            
            w.Read(w_in);
            
            if( Abs( w.Norm() - static_cast<Real>(1)) <= 128 * Scalar::eps<Real> )
            {
                // phi in [0, 1).
                const Real t = std::fmod(
                    std::atan2( w[1], w[0] ) * Scalar::TwoPiInv<Real> + static_cast<Real>(1),
                    static_cast<Real>(1)
                );
                
                // T in [0, edge_count_).
                const Real T = t * vertex_count_;
                const Real lambda = std::fmod(T,static_cast<Real>(1));
                const Int  j      = static_cast<Int>(std::floor(T));
                const Int  j_next = (j+1 < vertex_count_) ? j + 1 : 0;
                
                Vector_T a ( curve_, j      );
                Vector_T b ( curve_, j_next );
                
                Vector_T x; // Point on the unit sphere;
                
                LinearCombine( static_cast<Real>(1) - lambda, a, lambda, b, x );
                
                x.Normalize();
                
                x.Write(w_out);
                
                return;
            }
            
            w *= -static_cast<Real>(1);
            
            S2D_.ReadShiftVector( w.data() );
            
            S2D_.Shift();
            
            for( Int k = 0; k < edge_count_; ++k)
            {
                // Compute corresponding point on S^AmbDim by piecewise-linear interpolation.
                
                Vector2D_T y2D_k; // Point on the unit circle
                
                S2D_.WriteEdgeVectors(y2D_k,k);
                
                // phi in [0, 1).
                const Real t = std::fmod(
                    std::atan2( y2D_k[1], y2D_k[0] ) * Scalar::TwoPiInv<Real> + static_cast<Real>(1),
                    static_cast<Real>(1)
                );
                
                // T in [0, edge_count_).
                const Real T = t * vertex_count_;
                const Real lambda = std::fmod(T,static_cast<Real>(1));
                const Int  j      = static_cast<Int>(std::floor(T));
                const Int  j_next = (j+1 < vertex_count_) ? j + 1 : 0;
                
                Vector_T a ( curve_, j      );
                Vector_T b ( curve_, j_next );
                
                Vector_T x_k; // Point on the unit sphere;
                
                LinearCombine( static_cast<Real>(1) - lambda, a, lambda, b, x_k );
                
                S_.ReadInitialEdgeVectors(x_k,k,true);
            }
            
            S_.ComputeConformalClosure();
            
            S_.WriteShiftVector( w_out );
        }
        
        void operator()( 
            cptr<Real> w_in,
            const Int point_count,
            mptr<Real> w_out,
            const Int thread_count
        )
        {
            ParallelDo(
                [&,this]( const Int thread )
                {
                    DouadyEarleExtension<AmbDim,Real,Int> E ( *this );
                    
                    const Int i_begin = JobPointer( point_count, thread_count, thread     );
                    const Int i_end   = JobPointer( point_count, thread_count, thread + 1 );
                    
                    for( Int i = i_begin; i < i_end; ++i )
                    {
                        E( &w_in[2*i], &w_out[AmbDim*i] );
                    }
                },
                thread_count
            );
        }
        
                
    public:
        
        std::string ClassName() const
        {
            return std::string("CoBarS::DouadyEarleExtension") + "<" + ToString(AmbDim) + "," + TypeName<Real> + "," + TypeName<Int>  + "," + ">";
        }
        
    }; // class DouadyEarleExtension
    
} // namespace CoBarS

