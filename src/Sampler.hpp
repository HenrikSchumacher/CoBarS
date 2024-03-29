#pragma once

namespace CoBarS
{
    
    template<
        int AmbDim_,
        typename Real_    = double,
        typename Int_     = int_fast32_t,
        typename PRNG_T_  = Xoshiro256Plus,
        bool vectorizeQ   = true,
        bool zerofyfirstQ = false
    >
    class Sampler : public SamplerBase<AmbDim_,Real_,Int_>
    {
        ASSERT_FLOAT(Real_);
        ASSERT_REAL(Real_);
        ASSERT_INT(Int_);

    private:
        
        using Base_T = SamplerBase<AmbDim_,Real_,Int_>;

    public:

        using Real   = Real_;
        using Int    = Int_;
        using PRNG_T = PRNG_T_;
        
        static constexpr Int AmbDim = int_cast<Int>(AmbDim_);
        
        using Vector_T           = typename Base_T::Vector_T;
        using SquareMatrix_T     = typename Base_T::SquareMatrix_T;
        using SymmetricMatrix_T  = typename Base_T::SymmetricMatrix_T;

        using RandomVariable_T   = typename Base_T::RandomVariable_T;
        using RandomVariable_Ptr = std::shared_ptr<RandomVariable_T>;

        using Weights_T          = typename Base_T::Weights_T;
        using Setting_T          = typename Base_T::Setting_T;

        using VectorList_T       = Tiny::VectorList<AmbDim,Real,Int>;
        using Matrix_T           = Tensor2<Real,Int>;
    
    public:
        
        Sampler() = default;
        
        ~Sampler() = default;
        
        explicit Sampler(
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   Base_T( settings_ )
        ,   edge_count(edge_count_)
        ,   r   ( edge_count, one / edge_count )
        ,   rho ( edge_count, one )
        ,   total_r_inv ( one )
        {
            ComputeEdgeSpaceSamplingHelper();
            ComputeEdgeQuotientSpaceSamplingHelper();
            
            if constexpr ( vectorizeQ )
            {
                x = VectorList_T( edge_count     );
                y = VectorList_T( edge_count     );
                p = VectorList_T( edge_count + 1 );
            }
            else
            {
                x = Matrix_T( edge_count    , AmbDim );
                y = Matrix_T( edge_count    , AmbDim );
                p = Matrix_T( edge_count + 1, AmbDim );
            }
        }
        
        explicit Sampler(
            cptr<Real> r_in,
            cptr<Real> rho_in,
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   Base_T( settings_ )
        ,   edge_count(edge_count_)
        ,   r   ( edge_count )
        ,   rho ( rho_in, edge_count )
        {
            ComputeEdgeSpaceSamplingHelper();
            ComputeEdgeQuotientSpaceSamplingHelper();
            
            if constexpr ( vectorizeQ )
            {
                x = VectorList_T( edge_count     );
                y = VectorList_T( edge_count     );
                p = VectorList_T( edge_count + 1 );
            }
            else
            {
                x = Matrix_T( edge_count    , AmbDim );
                y = Matrix_T( edge_count    , AmbDim );
                p = Matrix_T( edge_count + 1, AmbDim );
            }
            
            ReadEdgeLengths(r_in);
        }
        
        
        // Copy constructor
        Sampler( const Sampler & other )
        :   Base_T( other )
        ,   edge_count( other.edge_count )
        ,   x(other.x)
        ,   y(other.y)
        ,   p(other.p)
        ,   r(other.r)
        ,   rho(other.rho)
        ,   total_r_inv(other.total_r_inv)
        ,   w(other.w)
        ,   F(other.F)
        ,   DF(other.DF)
        ,   L(other.L)
        ,   u(other.u)
        ,   z(other.z)
        ,   iter(other.iter)
        ,   squared_residual(other.squared_residual)
        ,   residual(other.residual)
        ,   edge_space_sampling_helper              ( other.edge_space_sampling_helper              )
        ,   edge_quotient_space_sampling_helper     ( other.edge_quotient_space_sampling_helper     )
        ,   edge_space_sampling_weight              ( other.edge_space_sampling_weight              )
        ,   edge_quotient_space_sampling_correction ( other.edge_quotient_space_sampling_correction )
        ,   lambda_min(other.lambda_min)
        ,   q(other.q)
        ,   errorestimator(other.errorestimator)
        ,   linesearchQ(other.linesearchQ)
        ,   succeededQ(other.succeededQ)
        ,   continueQ(other.continueQ)
        ,   ArmijoQ(other.ArmijoQ)
        ,   moments(other.moments)
        {
            LoadRandomVariables( other.F_list );
        }
        
        friend void swap( Sampler & A, Sampler & B ) noexcept
        {
            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
            using std::swap;

//            swap( static_cast<Base_T&>(A), static_cast<Base_T&>(B) );
            
            swap(A.edge_count,B.edge_count);
            swap(A.x,B.x);
            swap(A.y,B.y);
            swap(A.p,B.p);
            swap(A.r,B.r);
            swap(A.rho,B.rho);
            swap(A.total_r_inv,B.total_r_inv);
            swap(A.w,B.w);
            swap(A.F,B.F);
            swap(A.DF,B.DF);
            swap(A.L,B.L);
            swap(A.u,B.u);
            swap(A.z,B.z);

            swap(A.iter,             B.iter             );
            swap(A.squared_residual, B.squared_residual );
            swap(A.residual,         B.residual         );
            
            swap(A.edge_space_sampling_helper,                  B.edge_space_sampling_helper              );
            swap(A.edge_quotient_space_sampling_helper,         B.edge_quotient_space_sampling_helper     );
            swap(A.edge_space_sampling_weight,                  B.edge_space_sampling_weight              );
            swap(A.edge_quotient_space_sampling_correction,     B.edge_quotient_space_sampling_correction );
            
            swap(A.lambda_min,B.lambda_min);
            swap(A.q,B.q);
            swap(A.errorestimator,B.errorestimator);
            swap(A.linesearchQ,B.linesearchQ);
            swap(A.succeededQ,B.succeededQ);
            swap(A.continueQ,B.continueQ);
            swap(A.ArmijoQ,B.ArmijoQ);
            
            swap(A.moments,B.moments);
            swap(A.F_list,B.F_list);
        }
        
        // Copy assignment operator
        Sampler & operator=(Sampler other)
        {
            // copy-and-swap idiom
            // see https://stackoverflow.com/a/3279550/8248900 for details
            swap(*this, other);

            return *this;
        }

        /* Move constructor */
        Sampler( Sampler && other ) noexcept
        :   Base_T()
        {
            swap(*this, other);
        }
        
    protected:
        
        Int edge_count = 0;
        
        mutable PRNG_T random_engine;
        
        mutable std::normal_distribution<Real> normal_dist {zero,one};
        
        mutable std::uniform_real_distribution<Real> phi_dist {Scalar::Zero<Real>,Scalar::TwoPi<Real>};
        
        Setting_T settings;

        std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> x;
        std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> y;
        std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> p;
        
        
        Weights_T      r {0};
        Weights_T    rho {0};
        
        Real total_r_inv = one;
        
        Vector_T w;           // current point in hyperbolic space.
        Vector_T F;           // right hand side of Newton iteration.
        SymmetricMatrix_T DF; // nabla F(0) with respect to measure ys
        SymmetricMatrix_T L;  // storing Cholesky factor.
        Vector_T u;           // update direction
        
        Vector_T z;           // Multiple purpose buffer.
        
        Int iter = 0;
        
        Real squared_residual = 1;
        Real         residual = 1;

        Real edge_space_sampling_helper              = 1;
        Real edge_quotient_space_sampling_helper     = 1;
        
        Real edge_space_sampling_weight              = 0;
        Real edge_quotient_space_sampling_correction = 0;
        
        Real lambda_min = eps;
        Real q = one;
        Real errorestimator = infty;
        
        bool linesearchQ = true;    // Toggle line search.
        bool succeededQ  = false;   // Whether algorithm has succeded.
        bool continueQ   = true;    // Whether to continue with the main loop.
        bool ArmijoQ     = false;   // Whether Armijo condition was met last time we checked.
        
        static constexpr Real zero              = 0;
        static constexpr Real half              = 0.5;
        static constexpr Real one               = 1;
        static constexpr Real two               = 2;
        static constexpr Real three             = 3;
        static constexpr Real four              = 4;
        static constexpr Real eps               = std::numeric_limits<Real>::min();
        static constexpr Real infty             = std::numeric_limits<Real>::max();
        static constexpr Real small_one         = 1 - 16 * eps;
        static constexpr Real big_one           = 1 + 16 * eps;
        static constexpr Real g_factor          = 4;
        static constexpr Real g_factor_inv      = one/g_factor;
        static constexpr Real norm_threshold    = static_cast<Real>(0.99 * 0.99 + 16 * eps);
        static constexpr Real two_pi            = Scalar::TwoPi<Real>;

        
    private:
        
        // TODO: Add copy behavior to theses?
        
        mutable Tensor2<Real,Int> moments;
        mutable std::vector<std::shared_ptr<RandomVariable_T>> F_list;
        
    protected:
        
#include "Sampler/Optimization.hpp"

#include "Sampler/Reweighting.hpp"
        
#include "Sampler/Methods.hpp"

#include "Sampler/RandomOpenPolygons.hpp"
        
#include "Sampler/RandomClosedPolygons.hpp"
        
#include "Sampler/Sample.hpp"
        
#include "Sampler/ConformalClosures.hpp"
        
#include "Sampler/BinnedSample.hpp"
        
#include "Sampler/ConfidenceSample.hpp"
        
    public:
        
        // This routine seems to be somewhat slow; better not use it if you can avoid it.
        virtual Real InitialEdgeCoordinates( const Int k, const Int i ) const override
        {
            return (vectorizeQ ? x[i][k] : x[k][i]);
        }
        
        virtual Vector_T InitialEdgeCoordinates( const Int k ) const override
        {
            return Vector_T ( x, k );
        }
        
        virtual void ReadInitialEdgeCoordinates( cptr<Real> x_in, bool normalize = true ) override
        {
            if( normalize )
            {
                for( Int k = 0; k < edge_count; ++k )
                {
                    Vector_T x_k ( x_in, k );
                    
                    x_k.Normalize();
                    
                    x_k.Write( x, k );
                }
            }
            else
            {
                x.Read( x_in );
            }
        }

        void ReadInitialEdgeCoordinates( cref<Vector_T> x_in, const Int k, bool normalize = true )
        {
            Vector_T x_k ( x_in );
            
            if( normalize )
            {
                x_k.Normalize();
            }
            
            x_k.Write( x, k );
        }
        
        virtual void ReadInitialEdgeCoordinates( cptr<Real> x_in, const Int k, bool normalize = true ) override
        {
            ReadInitialEdgeCoordinates( &x_in[ AmbDim * edge_count * k], normalize );
        }
        
        virtual void WriteInitialEdgeCoordinates( mptr<Real> x_out ) const override
        {
            x.Write( x_out );
        }
        
        virtual void WriteInitialEdgeCoordinates( mptr<Real> x_out, const Int k ) const override
        {
            WriteInitialEdgeCoordinates( &x_out[ AmbDim * edge_count * k ]);
        }

        void DumpInitialEdgeCoordinates() const
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                valprint( "x["+ToString(i)+"]", x[i] );
            }
        }
        
        void DumpEdgeCoordinates() const
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                valprint( "y["+ToString(i)+"]", y[i] );
            }
        }
        
        virtual void RandomizeInitialEdgeCoordinates() override
        {
            for( Int k = 0; k < edge_count; ++k )
            {
                Vector_T x_k;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    x_k[i] = normal_dist( random_engine );
                }
                
                x_k.Normalize();
                
                x_k.Write( x, k );
            }
        }
        
        
        virtual void ComputeConformalClosure() override
        {
            ComputeShiftVector();

            Optimize();
        }
        
        // This routine seems to be somewhat slow; better not use it if you can avoid it.
        virtual Real EdgeCoordinates( const Int k, const Int i ) const override
        {
            return COND( vectorizeQ, y[i][k], y[k][i] );
        }
        
        virtual Vector_T EdgeCoordinates( const Int k ) const override
        {
            return Vector_T ( y, k );
        }
        
        virtual void ReadEdgeCoordinates( cptr<Real> y_in ) override
        {
            y.Read( y_in );
        }
        
        virtual void ReadEdgeCoordinates( cptr<Real> y_in, const Int k ) override
        {
            ReadEdgeCoordinates( &y_in[ AmbDim * edge_count * k ]);
        }
        
        virtual void WriteEdgeCoordinates( mptr<Real> y_out ) const override
        {
            y.Write( y_out );
        }
        
        virtual void WriteEdgeCoordinates( mptr<Real> y_out, const Int k ) const override
        {
            WriteEdgeCoordinates( &y_out[ AmbDim * edge_count * k ]);
        }
        
        void WriteEdgeCoordinates( mref<Vector_T> y_out, const Int k ) const
        {
            y_out.Read( y, k );
        }
        

        // This routine seems to be somewhat slow; better not use it if you can avoid it.
        virtual Real SpaceCoordinates( const Int k, const Int i ) const override
        {
            return COND(vectorizeQ, p[i][k], p[k][i] );
        }
        
        virtual Vector_T SpaceCoordinates( const Int k ) const override
        {
            return Vector_T ( p, k );
        }
        
        virtual void WriteSpaceCoordinates( mptr<Real> p_out ) const override
        {
            p.Write( p_out );
        }
        
        virtual void WriteSpaceCoordinates( mptr<Real> p_out, const Int k ) const override
        {
            WriteSpaceCoordinates( &p_out[ (edge_count+1) * AmbDim * k ]);
        }
        
        virtual void ComputeSpaceCoordinates() override
        {
            //Caution: This gives only half the weight to the end vertices of the chain.
            //Thus this is only really the barycenter, if the chain is closed!
            
            // We treat the edges as massless.
            // All mass is concentrated in the vertices, and each vertex carries the same mass.
            
            Vector_T barycenter        (zero);
            Vector_T point_accumulator (zero);
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Vector_T y_k ( y, k );
                
                const Real r_k = r[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real offset = r_k * y_k[i];
                    
                    barycenter[i] += (point_accumulator[i] + half * offset);
                    
                    point_accumulator[i] += offset;
                }
            }
            
            barycenter *= Inv<Real>( edge_count );
            
            point_accumulator = barycenter;
            
            point_accumulator *= -one;
            
            point_accumulator.Write( p, 0 );
  
            for( Int k = 0; k < edge_count; ++k )
            {
                const Vector_T y_k ( y, k );
                
                const Real r_k = r[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    point_accumulator[i] += r_k * y_k[i];
                }
                
                point_accumulator.Write( p, k + 1 );
            }
        }
        
        virtual const std::vector<std::shared_ptr<RandomVariable_T>> & RandomVariables() const override
        {
            return F_list;
        }
        
        virtual Int RandomVariablesCount() const override
        {
            return static_cast<Int>(F_list.size());
        }
        
        virtual void LoadRandomVariables( const std::vector<std::shared_ptr<RandomVariable_T>> & F_list_ ) const override
        {
            F_list.clear();
            
            for( RandomVariable_Ptr f : F_list_ )
            {
                F_list.push_back( f->Clone() );
            }
        }
        
        virtual void ClearRandomVariables() const override
        {
            F_list.clear();
        }

        virtual Real EvaluateRandomVariable( Int i ) const override
        {
            return (*F_list[i])( *this );
        }

        
    public:
        
        
        std::string PRNG_Name() const override
        {
//            return random_engine[0].ClassName();
            return random_engine.ClassName();
        }
        
        std::string ClassName() const override
        {
            return std::string("CoBarS::Sampler") + "<" + ToString(AmbDim) + "," + TypeName<Real> + "," + TypeName<Int>  + "," + PRNG_Name() + "," + ToString(vectorizeQ) +"," + ToString(zerofyfirstQ) + ">";
        }
        
    }; // class Sampler
    
} // namespace CoBarS
