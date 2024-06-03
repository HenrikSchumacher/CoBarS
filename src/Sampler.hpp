#pragma once

namespace CoBarS
{
    
/*!
 * @brief The main class of CoBarS. It does the actual sampling.
 *
 * @tparam AmbDim_ The dimension of the ambient space.
 *
 * @tparam Real_ A real floating point type.
 *
 * @tparam Int_  An integer type.
 *
 * @tparam PRNG_T_ A class of a pseudorandom number generator.
 * Possible values are
 *  - CoBarS::MT64
 *  - CoBarS::PCG64
 *  - CoBarS::WyRand
 *  - CoBarS::Xoshiro256Plus
 */
    
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
        static_assert(FloatQ<Real_>,"");
        static_assert(Scalar::RealQ<Real_>,"");
        static_assert(IntQ<Int_>,"");

    private:
        
        using Base_T = SamplerBase<AmbDim_,Real_,Int_>;

    public:
        
        using typename Base_T::Real;
        using typename Base_T::Int;
        using PRNG_T = PRNG_T_;
        
        using Base_T::AmbDim;
        
        using typename Base_T::Vector_T;
        using typename Base_T::SquareMatrix_T;
        using typename Base_T::SymmetricMatrix_T;

        using typename Base_T::Weights_T;
        using typename Base_T::Setting_T;
        
        
        using typename Base_T::RandomVariable_T;
        using RandomVariable_Ptr = std::shared_ptr<RandomVariable_T>;

        using VectorList_T = Tiny::VectorList<AmbDim,Real,Int>;
        using Matrix_T     = Tensor2<Real,Int>;
    
    public:
        
        Sampler() = default;
        
        ~Sampler() = default;
        
        explicit Sampler(
            const Int edge_count,
            const Setting_T settings = Setting_T()
        )
        :   Base_T( settings )
        ,   edge_count_(edge_count)
        ,   r_   ( edge_count_, one / edge_count_ )
        ,   rho_ ( edge_count_, one )
        ,   total_r_inv ( one )
        {
            ComputeEdgeSpaceSamplingHelper();
            ComputeEdgeQuotientSpaceSamplingHelper();
            
            if constexpr ( vectorizeQ )
            {
                x_ = VectorList_T( edge_count_     );
                y_ = VectorList_T( edge_count_     );
                p_ = VectorList_T( edge_count_ + 1 );
            }
            else
            {
                x_ = Matrix_T( edge_count_    , AmbDim );
                y_ = Matrix_T( edge_count_    , AmbDim );
                p_ = Matrix_T( edge_count_ + 1, AmbDim );
            }
        }
        
        explicit Sampler(
            const Real * restrict const r,
            const Real * restrict const rho,
            const Int edge_count,
            const Setting_T settings = Setting_T()
        )
        :   Base_T      ( settings )
        ,   edge_count_ (edge_count)
        ,   r_          ( edge_count_ )
        ,   rho_        ( rho, edge_count_ )
        {
            ComputeEdgeSpaceSamplingHelper();
            ComputeEdgeQuotientSpaceSamplingHelper();
            
            if constexpr ( vectorizeQ )
            {
                x_ = VectorList_T( edge_count_     );
                y_ = VectorList_T( edge_count_     );
                p_ = VectorList_T( edge_count_ + 1 );
            }
            else
            {
                x_ = Matrix_T( edge_count_    , AmbDim );
                y_ = Matrix_T( edge_count_    , AmbDim );
                p_ = Matrix_T( edge_count_ + 1, AmbDim );
            }
            
            ReadEdgeLengths(r);
        }
        
        /*!
         * @brief Copy constructor.
         */
        
        Sampler( const Sampler & other )
        :   Base_T( other )
        ,   edge_count_( other.edge_count_ )
        ,   x_(other.x_)
        ,   y_(other.y_)
        ,   p_(other.p_)
        ,   p_initializedQ(other.p_initializedQ)
        ,   r_(other.r_)
        ,   rho_(other.rho_)
        ,   total_r_inv(other.total_r_inv)
        ,   w_(other.w_)
        ,   F_(other.F_)
        ,   DF_(other.DF_)
        ,   L(other.L)
        ,   u_(other.u_)
        ,   z_(other.z_)
        ,   iter(other.iter)
        ,   squared_residual(other.squared_residual)
        ,   residual(other.residual)
        ,   edge_space_sampling_helper              ( other.edge_space_sampling_helper              )
        ,   edge_quotient_space_sampling_helper     ( other.edge_quotient_space_sampling_helper     )
        ,   edge_space_sampling_weight              ( other.edge_space_sampling_weight              )
        ,   edge_quotient_space_sampling_weight     ( other.edge_quotient_space_sampling_weight )
        ,   lambda_min(other.lambda_min)
        ,   q(other.q)
        ,   errorestimator(other.errorestimator)
        ,   linesearchQ(other.linesearchQ)
        ,   succeededQ(other.succeededQ)
        ,   continueQ(other.continueQ)
        ,   ArmijoQ(other.ArmijoQ)
        ,   moments_(other.moments_)
        {
            LoadRandomVariables( other.F_list_ );
        }
        
        /*!
         * @brief Swap routine; see [stackoverflow.com](https://stackoverflow.com/questions/5695548/public-friend-swap-member-function) for details.
         */

        friend void swap( Sampler & A, Sampler & B ) noexcept
        {
            using std::swap;
   
            swap(A.edge_count_,B.edge_count_);
            swap(A.x_,B.x_);
            swap(A.y_,B.y_);
            swap(A.p_,B.p_);
            swap(A.p_initializedQ,B.p_initializedQ);
            swap(A.r_,B.r_);
            swap(A.rho_,B.rho_);
            swap(A.total_r_inv,B.total_r_inv);
            swap(A.w_,B.w_);
            swap(A.F_,B.F_);
            swap(A.DF_,B.DF_);
            swap(A.L,B.L);
            swap(A.u_,B.u_);
            swap(A.z_,B.z_);

            swap(A.iter,             B.iter             );
            swap(A.squared_residual, B.squared_residual );
            swap(A.residual,         B.residual         );
            
            swap(A.edge_space_sampling_helper,          B.edge_space_sampling_helper          );
            swap(A.edge_quotient_space_sampling_helper, B.edge_quotient_space_sampling_helper );
            swap(A.edge_space_sampling_weight,          B.edge_space_sampling_weight          );
            swap(A.edge_quotient_space_sampling_weight, B.edge_quotient_space_sampling_weight );
            
            swap(A.lambda_min,B.lambda_min);
            swap(A.q,B.q);
            swap(A.errorestimator,B.errorestimator);
            swap(A.linesearchQ,B.linesearchQ);
            swap(A.succeededQ,B.succeededQ);
            swap(A.continueQ,B.continueQ);
            swap(A.ArmijoQ,B.ArmijoQ);
            
            swap(A.moments_,B.moments_);
            swap(A.F_list_,B.F_list_);
        }
        
        /*!
         * @brief Copy assignment operator, using [copy-and-swap idiom](https://stackoverflow.com/a/3279550/8248900).
         */
        
        Sampler & operator=(Sampler other)
        {
            swap(*this, other);

            return *this;
        }

        /*!
         * @brief Move constructor, using [copy-and-swap idiom](https://stackoverflow.com/a/3279550/8248900).
         */
        
        /* Move constructor */
        Sampler( Sampler && other ) noexcept
        :   Base_T()
        {
            swap(*this, other);
        }
        
    private:
        
        /*!
         * @brief Number of edges in the represented polygon.
         */
        
        Int edge_count_ = 0;
        
        /*!
         * @brief The instance's own pseudorandom number generator.
         */
        
        mutable PRNG_T random_engine;
        
        mutable std::normal_distribution<Real> normal_dist {zero,one};
        
        mutable std::uniform_real_distribution<Real> phi_dist {0,Scalar::TwoPi<Real>};
        
        using Base_T::settings_;

        /*!
         * @brief The open polyline's unit edge vectors.
         */
        
        std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> x_;
        
        /*!
         * @brief The shifted polyline's unit edge vectors. (After the optimization has succeeded: the closed polygon's unit edge vectors.)
         */
        
        std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> y_;
        
        /*!
         * @brief The closed polyline's vertex coordinates.
         */
        mutable std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> p_;
        
        /*!
         * @brief Boolean that indicates whether `p_initializedQ` has been recomputed already.
         */
        mutable bool p_initializedQ = false;
        
        
        /*!
         * @brief The edge lengths of the represented polygon.
         */
        
        Weights_T r_ {0};
        
        /*!
         * @brief The weights of the Riemannian metric used on the Cartesian product of unit spheres.
         */
        
        Weights_T rho_ {0};
        
        /*!
         * @brief Inverse of the total arc length of the represented polygon.
         */
        
        Real total_r_inv = one;
        
        /*!
         * @brief The current shift vector in hyperbolic space. (After the optimization has succeeded: the conformal barycenter.)
         */
        
        Vector_T w_;
        
        /*!
         * @brief Right hand side of Newton iteration.
         */
        
        Vector_T F_;
        
        /*!
         * @brief Nabla F(0) with respect to probability measure defined by point cloud `y_` and weights `r_`.
         */
        
        SymmetricMatrix_T DF_;
        
        /*!
         * @brief An auxiliary buffer to store Cholesky factors.
         */
        SymmetricMatrix_T L;
        
        /*!
         * @brief Update direction.
         */
        Vector_T u_;
        
        
        /*!
         * @brief Multiple purpose buffer.
         */
        Vector_T z_;
        
        /*!
         * @brief Number of optimization iterations used so far.
         */
        
        Int iter = 0;
        
        Real squared_residual = 1;
        Real         residual = 1;

        Real edge_space_sampling_helper                  = 1;
        Real edge_quotient_space_sampling_helper         = 1;
        
        mutable Real edge_space_sampling_weight          = -1;
        mutable Real edge_quotient_space_sampling_weight = -1;
        
        Real lambda_min = eps;
        Real q = one;
        Real errorestimator = infty;
        
        bool linesearchQ = true;    // Toggle line search.
        bool succeededQ  = false;   // Whether algorithm has succeded.
        bool continueQ   = true;    // Whether to continue with the main loop.
        bool ArmijoQ     = false;   // Whether Armijo condition was met last time we checked.
        
    private:
        
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
        static constexpr Real norm_threshold    = 0.99 * 0.99 + 16 * eps;
        static constexpr Real two_pi            = Scalar::TwoPi<Real>;

        
        // TODO: Add copy behavior to these?
        
        mutable Tensor2<Real,Int> moments_;
        mutable std::vector<std::shared_ptr<RandomVariable_T>> F_list_;
        
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
        
//        virtual Real InitialEdgeCoordinate( const Int k, const Int i ) const override
//        {
//            return (vectorizeQ ? x_[i][k] : x_[k][i]);
//        }
        
        virtual Vector_T InitialEdgeVector( const Int k ) const override
        {
            return Vector_T (x_,k);
        }
        
//        virtual void ReadInitialEdgeVector(
//            const Vector_T & x, const Int k, bool normalizeQ = true
//        ) override
//        {
//            Vector_T x_k (x);
//            
//            if( normalizeQ )
//            {
//                x_k.Normalize();
//            }
//            
//            x_k.Write(x_,k);
//        }

        virtual void ReadInitialEdgeVectors(
            const Real * restrict const x, const Int offset = 0 , bool normalizeQ = true
        ) override
        {
            cptr<Real> X = &x[AmbDim * edge_count_ * offset];
            
            if( normalizeQ )
            {
                for( Int k = 0; k < edge_count_; ++k )
                {
                    Vector_T x_k (X,k);
                    
                    x_k.Normalize();
                    
                    x_k.Write(x_,k);
                }
            }
            else
            {
                x_.Read(X);
            }
        }
        
        virtual void WriteInitialEdgeVectors( 
            Real * restrict const x, const Int offset = 0
        ) const override
        {
            x_.Write( &x[ AmbDim * edge_count_ * offset ]);
        }

        void DumpInitialEdgeVectors() const
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                valprint( "x_["+ToString(i)+"]", x_[i] );
            }
        }
        
        void DumpEdgeVectors() const
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                valprint( "y_["+ToString(i)+"]", y_[i] );
            }
        }
        
        virtual void RandomizeInitialEdgeVectors() override
        {
            for( Int k = 0; k < edge_count_; ++k )
            {
                Vector_T x_k;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    x_k[i] = normal_dist( random_engine );
                }
                
                x_k.Normalize();
                
                x_k.Write(x_,k);
            }
        }
        
        virtual void ComputeConformalClosure() override
        {
            ResetResults();
            
            ComputeInitialShiftVector();

            Optimize();
        }
        

        
        virtual Vector_T EdgeVector( const Int k ) const override
        {
            return Vector_T (y_,k);
        }
        
        virtual void ReadEdgeVectors( const Real * restrict const y, const Int offset = 0 ) override
        {
            ResetResults();
            
            y_.Read( &y[ AmbDim * edge_count_ * offset ]);
        }
        
        virtual void WriteEdgeVectors( Real * restrict const y, const Int offset = 0 ) const override
        {
            y_.Write( &y[ AmbDim * edge_count_ * offset ]);
        }
        

//        // This routine seems to be somewhat slow; better not use it if you can avoid it.
//        virtual Real VertexCoordinate( const Int k, const Int i ) const override
//        {
//            RequireVertexPositions();
//
//            return COND(vectorizeQ, p_[i][k], p_[k][i] );
//        }
        
        virtual Vector_T VertexPosition( const Int k ) const override
        {
            RequireVertexPositions();

            return Vector_T (p_,k);
        }
        
        virtual void WriteVertexPositions( Real * restrict const p, const Int offset = 0 ) const override
        {
            RequireVertexPositions();
            
            p_.Write( &p[ (edge_count_+1) * AmbDim * offset ] );
        }
        
        virtual void RequireVertexPositions() const override
        {
            //Caution: This gives only half the weight to the end vertices of the chain.
            //Thus this is only really the barycenter, if the chain is closed!
            
            // We treat the edges as massless.
            // All mass is concentrated in the vertices, and each vertex carries the same mass.
            
            if( !p_initializedQ )
            {
                Vector_T barycenter        (zero);
                Vector_T point_accumulator (zero);
                
                for( Int k = 0; k < edge_count_; ++k )
                {
                    const Vector_T y_k ( y_, k );
                    
                    const Real r_k = r_[k];
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        const Real offset = r_k * y_k[i];
                        
                        barycenter[i] += (point_accumulator[i] + half * offset);
                        
                        point_accumulator[i] += offset;
                    }
                }
                
                barycenter *= Inv<Real>( edge_count_ );
                
                point_accumulator = barycenter;
                
                point_accumulator *= -one;
                
                point_accumulator.Write( p_, 0 );
                
                for( Int k = 0; k < edge_count_; ++k )
                {
                    const Vector_T y_k ( y_, k );
                    
                    const Real r_k = r_[k];
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        point_accumulator[i] += r_k * y_k[i];
                    }
                    
                    point_accumulator.Write( p_, k + 1 );
                }
            }
        }
        
        virtual const std::vector<std::shared_ptr<RandomVariable_T>> & RandomVariables() const override
        {
            return F_list_;
        }
        
        virtual Int RandomVariablesCount() const override
        {
            return static_cast<Int>(F_list_.size());
        }
        
        virtual void LoadRandomVariables( const std::vector<std::shared_ptr<RandomVariable_T>> & F_list ) const override
        {
            F_list_.clear();
            
            for( RandomVariable_Ptr F : F_list )
            {
                F_list_.push_back( F->Clone() );
            }
        }
        
        virtual void ClearRandomVariables() const override
        {
            F_list_.clear();
        }

        virtual Real EvaluateRandomVariable( Int i ) const override
        {
            return (*F_list_[i])( *this );
        }

        
    public:
        
        
        std::string PRNG_Name() const override
        {
            return random_engine.ClassName();
        }
        
        std::string ClassName() const override
        {
            return std::string("CoBarS::Sampler") + "<" + ToString(AmbDim) + "," + TypeName<Real> + "," + TypeName<Int>  + "," + PRNG_Name() + "," + ToString(vectorizeQ) +"," + ToString(zerofyfirstQ) + ">";
        }
        
    }; // class Sampler
    
} // namespace CoBarS
