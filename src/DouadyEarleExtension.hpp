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
        
        ASSERT_FLOAT(Real_);
        ASSERT_REAL(Real_);
        ASSERT_INT(Int_);
        
    public:

        using Real   = Real_;
        using Int    = Int_;
        
        static constexpr Int AmbDim = int_cast<Int>(AmbDim_);
        
        // It does not matter which pseudorandom number generator we use.
        // Nothing is random here.
        using PRNG_T  = Xoshiro256Plus;
        
        using Sampler2D_T        = Sampler<2     ,Real,Int,PRNG_T,true,true>;
        using Sampler_T          = Sampler<AmbDim,Real,Int,PRNG_T,true,true>;

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
        

        const Int edge_count;
        
        mutable Sampler2D_T S2D;
        mutable Sampler_T   S;
        
        mutable Int vertex_count;
        mutable VectorList_T curve;
        
    public:
        
        DouadyEarleExtension() = default;
        
        ~DouadyEarleExtension() = default;
        
        explicit DouadyEarleExtension(
            const Int edge_count_
        )
        :   edge_count ( edge_count_  )
        ,   S2D        ( edge_count )
        ,   S          ( edge_count )
        {
            const Real step  = Scalar::TwoPi<Real> / edge_count;
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Real angle = step * k;
                
                Vector2D_T z = {std::cos(angle), std::sin(angle) };
                
                S2D.ReadInitialEdgeCoordinates( z, k, false );
            }
        }
        
        
        // Copy constructor
        DouadyEarleExtension( const DouadyEarleExtension & other )
        :   edge_count   ( other.edge_count   )
        ,   S2D          ( other.S2D          )
        ,   S            ( other.S            )
        ,   vertex_count ( other.vertex_count )
        ,   curve        ( other.curve        )
        {}
        
        void LoadCurve( cptr<Real> curve_, const Int vertex_count_,
            const bool normalizeQ = false
        )
        {
            vertex_count = vertex_count_;
            
            if( curve.Dimension(1) != vertex_count )
            {
                curve = VectorList_T( vertex_count );
            }
            
            if( normalizeQ )
            {
                for( Int k = 0; k < vertex_count; ++k )
                {
                    Vector_T u;
                    u.Read( &curve_[AmbDim * k] );
                    u.Normalize();
                    u.Write( curve, k );
                }
            }
            else
            {
                for( Int k = 0; k < vertex_count; ++k )
                {
                    Vector_T u;
                    u.Read( &curve_[AmbDim * k] );
                    u.Write( curve, k );
                }
            }
        }
        
        void operator()( cptr<Real> w_in, mptr<Real> w_out, bool UseOldResultQ = false )
        {
            // Take 2D input vector and shift it to 0 so that
            Vector2D_T w;
            
            w.Read( w_in );
            
            if( Abs( w.Norm() - Scalar::One<Real>) <= 128 * Scalar::eps<Real> )
            {
                
                // phi in [0, 1).
                const Real t = std::fmod(
                    std::atan2( w[1], w[0] ) * Scalar::TwoPiInv<Real> + Scalar::One<Real>,
                    Scalar::One<Real>
                );
                
                // T in [0, edge_count).
                const Real T = t * vertex_count;
                const Real lambda = std::fmod(T,Scalar::One<Real>);
                const Int  j      = static_cast<Int>(std::floor(T));
                const Int  j_next = (j+1 < vertex_count) ? j + 1 : 0;
                
                Vector_T a ( curve, j      );
                Vector_T b ( curve, j_next );
                
                Vector_T x; // Point on the unit sphere;
                
                LinearCombine( Scalar::One<Real> - lambda, a, lambda, b, x );
                
                x.Normalize();
                
                x.Write( w_out );
                
                return;
            }
            
            w *= -Scalar::One<Real>;
            
            S2D.ReadShiftVector( w.data() );
            
            S2D.Shift();
            
            for( Int k = 0; k < edge_count; ++k)
            {
                // Compute corresponding point on S^AmbDim by piecewise-linear interpolation.
                
                Vector2D_T y2D_k; // Point on the unit circle
                
                S2D.WriteEdgeCoordinates( y2D_k, k );
                
//                if( k == 0 ) dump(y2D_k);
                
                // phi in [0, 1).
                const Real t = std::fmod(
                    std::atan2( y2D_k[1], y2D_k[0] ) * Scalar::TwoPiInv<Real> + Scalar::One<Real>,
                    Scalar::One<Real>
                );
                
                // T in [0, edge_count).
                const Real T = t * vertex_count;
                const Real lambda = std::fmod(T,Scalar::One<Real>);
                const Int  j      = static_cast<Int>(std::floor(T));
                const Int  j_next = (j+1 < vertex_count) ? j + 1 : 0;
                
                Vector_T a ( curve, j      );
                Vector_T b ( curve, j_next );
                
                Vector_T x_k; // Point on the unit sphere;
                
                LinearCombine( Scalar::One<Real> - lambda, a, lambda, b, x_k );
                
                S.ReadInitialEdgeCoordinates( x_k, k, true );
            }
            
//            S.DumpInitialEdgeCoordinates();
            
            if( UseOldResultQ )
            {
                // We use the shift vector that is still stored in S.
                // Should work well if w_in is close to the previous w_in.
            }
            else
            {
                S.ComputeShiftVector();
            }
            

            
            S.Optimize();
            
            S.WriteShiftVector( w_out );
        }
        
        void operator()( 
            cptr<Real> w_in,
            const Int point_count,
            mptr<Real> w_out,
            const Int thread_count,
            bool UseOldResultQ = false
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
                        E( &w_in[2*i], &w_out[AmbDim*i], UseOldResultQ);
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

