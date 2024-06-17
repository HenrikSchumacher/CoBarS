#pragma once

namespace AAM
{
    using namespace Tools;
    using namespace Tensors;
    
    /*!
     * @brief Implements the _(Progressive) Action-Angle Method_.
     *
     * @tparam REAL A real floating point type.
     *
     * @tparam INT  An integer type.
     *
     * @tparam PRNG_T A class of a pseudorandom number generator.
     * Possible values are
     *  - CoBarS::MT64
     *  - CoBarS::PCG64
     *  - CoBarS::WyRand
     *  - CoBarS::Xoshiro256Plus
     *
     *  @tparam PROGRESSIVE_Q If set to `false`, it uses the _Action-Angle Method_ by (Cantarella, Duplantier, Shonkwiler, and Uehara)[https://iopscience.iop.org/article/10.1088/1751-8113/49/27/275202]. Otherwise, it uses the the _Progressive Action-Angle Method_ sampler by Cantarella, Schumacher, and Shonkwiler.
     */
    
    template<
        typename REAL    = double,
        typename INT     = std::int_fast32_t,
        typename PRNG_T  = CoBarS::Xoshiro256Plus,
        bool PROGRESSIVE_Q = true
    >
    class Sampler
    {
        static_assert(FloatQ<REAL>,"");
        static_assert(IntQ<INT>,"");
        
    public:
        
        using Real   = REAL;
        using Int    = INT;
        using Prng_T = PRNG_T;
        
        using SpherePoints_T = Tensor2<Real,Int>;
        using SpacePoints_T  = Tensor2<Real,Int>;
        
        static constexpr Int AmbDim = 3;
        static constexpr bool ProgressiveQ = PROGRESSIVE_Q;
        
        using Vector_T = Tensors::Tiny::Vector<AmbDim,Real,Int>;
        
        Sampler() = default;
        
        ~Sampler(){}
        
        explicit Sampler( const Int edge_count )
        :   edge_count_ ( edge_count )
        {}

    protected:
        
        static constexpr Real zero              = 0;
        static constexpr Real half              = 0.5;
        static constexpr Real one               = 1;
        static constexpr Real two               = 2;
        static constexpr Real three             = 3;
        static constexpr Real four              = 4;
        static constexpr Real eps               = std::numeric_limits<Real>::min();
        static constexpr Real infty             = std::numeric_limits<Real>::max();
        static constexpr Real small_one         = one - 16 * eps;
        static constexpr Real big_one           = one + 16 * eps;
//        static constexpr Real two_pi            = Scalar::TwoPi<Real>;
        
        const Int edge_count_;
        
        mutable Prng_T random_engine;
                    
        std::uniform_real_distribution<Real> dist_1 { -one, one };
        
        std::uniform_real_distribution<Real> dist_2 { -Scalar::Pi<Real>, Scalar::Pi<Real> };
        
//        Tools::uniform_dist<Real,UInt64> dist_1 {-one, one};
//        
//        Tools::uniform_dist<Real,UInt64> dist_2 { -Scalar::Pi<Real>, Scalar::Pi<Real> };
        
        static constexpr Int max_trials = 10000000;
        
    public:
        
        /*!
         * Generates a single closed, equilateral random polygon.
         *
         * This is a port of the routine `plc_random_equilateral_closed_polygon` from the C library [plCurve](https://jasoncantarella.com/wordpress/software/plcurve) by Ted Ashton, Jason Cantarella, Harrison Chapman, and Tom Eddy.
         *
         * @param p Buffer for vertex positions; assumed to be of size `EdgeCount() * AmbientDimension()`. Coordinates are stored in interleaved form, i.e. the coordinates of each vertex lie contiguously in memory.
         */
        
        Int CreateRandomClosedPolygon( Real * restrict const p )
        {
            const Int n = edge_count_;
            
            // We use the user-supplied buffer as scratch space for the diagonal lengths.
            mptr<Real> d = &p[2*edge_count_+1];
            
            bool rejectedQ = true;
            
            Int trials = 0;
            
            // We need to find d[0], d[1], ... , d[n-2] from the moment polytope.
            // d[0], d[1], ..., d[n-2] are diagonals of a fan polyhedron if and only if the following is satisfied:
            //
            // (1)      d[0  ] = 1;
            // (2)      d[n-2] = 1;
            
            // (4)      |d[i-1] - d[i]| <= 1    for i = 1,...,n-3
            // (5)      d[i-1] + d[i] >= 1      for i = 1,...,n-3
            // (6)      d[i] > 0                for i = 1,...,n-3 -- obsolete, follows from (4) and (5)!
            
            // (7)      |d[n-3] - d[n-2]| <= 1
            // (8)      d[n-3] + d[n-2] >= 1
            
            // Because d[n-2] == 1, the latter two reduce to
            
            // (7')     |d[n-3] - 1| <= 1
            // (8')     d[n-3] >= 0             obsolete, follows from (6) for i = n-3.
            
            // Because d[n-3] >0, the condition (4') reduces to
            
            //  (7'') d[n-3] <= 2.
            
            // For condition (7'') it is necessary that d[i] - (n - 3 - i ) <= d[n-3] <= 2.
            // So we may stop early if
            //
            // d[i] - (n - 3 - i ) > 2.
            //
            // which is equivalent to
            //
            //  (7''')  d[i] + i > n - 1.
            //
            // However, this has a chance to be satisfied only in step i roughly (n - 2)/2 or later.
            // Anyways, checking this comes at an extra cost; so (7''') does not seem to help.
            
            // This guarantees (1), (2).
            d[0  ] = one;
            d[n-2] = one;
            
            while( rejectedQ && (trials < max_trials) )
            {
                ++trials;
                
                if constexpr ( ProgressiveQ )
                {
                    for( Int i = 1; i < (n-2); ++i )
                    {
                        // This guarantees (4):
                        d[i] = d[i-1] + dist_1( random_engine );
                        
                        // Check condition (5).
                        rejectedQ = d[i-1] + d[i] < one;
                        
                        if( rejectedQ )
                        {
                            break;
                        }
                    }
                }
                else
                {
                    for( Int i = 1; i < (n-2); ++i )
                    {
                        // This guarantees (4):
                        d[i] = d[i-1] + dist_1( random_engine );
                    }
                    
                    rejectedQ = false;
                    
                    // Check condition (5.
                    for( Int i = 1; i < (n-2); ++i )
                    {
                        rejectedQ = rejectedQ || (d[i-1] + d[i] < one);
                    }
                }
                
                // Check condition (7'') for last diagonal:
                rejectedQ = rejectedQ || ( d[n-3] > two );
            }
            
            if( trials > max_trials )
            {
                eprint("No success after " + ToString(max_trials) + " trials.");
                return trials;
            }
            
            p[0] = zero;
            p[1] = zero;
            p[2] = zero;
            p[3] = one ;
            p[4] = zero;
            p[5] = zero;
            
            // Unit vector pointing to position of currect vertex..
            Vector_T e { one, zero, zero };
            
            // Some buffer for the cross product in Rodrigues' formula.
            Vector_T v;
            
            // The current triangle's normal.
            Vector_T nu { zero, zero, one };
                        
            for( Int i = 0; i < n - 3; ++i )
            {
                // Next we compute the new unit vector that points to e by rotating e by the angle alpha about the unit normal of the triangle.
                
                // Compute angle alpha beteen diagonals d[i-1] and d[i] by the cosine identity.
                const Real cos_alpha = ( d[i] * d[i] + d[i+1] * d[i+1] - one )/( two * d[i] * d[i+1] );
                
                // 0 < alpha < Pi, so sin(alpha) is postive. Thus the following is safe.
                const Real sin_alpha = std::sqrt(one - cos_alpha * cos_alpha);
                
                Cross(nu,e,v);
                
                const Real factor = Dot(nu,e) * (one-cos_alpha);
                
                // Apply Rodrigues' formula
                e[0] = e[0] * cos_alpha + v[0] * sin_alpha + nu[0] * factor;
                e[1] = e[1] * cos_alpha + v[1] * sin_alpha + nu[1] * factor;
                e[2] = e[2] * cos_alpha + v[2] * sin_alpha + nu[2] * factor;
                
                // Normalize for stability
                e.Normalize();
                
                // Compute the new vertex position.
                p[3 * (i+2) + 0] = d[i+1] * e[0];
                p[3 * (i+2) + 1] = d[i+1] * e[1];
                p[3 * (i+2) + 2] = d[i+1] * e[2];
                
                // Now we also rotate the triangle's unit normal by theta[i] about the unit vector e.
                
                Cross(e,nu,v);
                
                const Real theta_i   = dist_2( random_engine );
                const Real cos_theta = std::cos(theta_i);
                const Real sin_theta = std::sin(theta_i);
                const Real factor_2  = Dot(e,nu) * (one-cos_theta);
                
                // Apply Rodrigues' formula
                nu[0] = nu[0] * cos_theta + v[0] * sin_theta + e[0] * factor_2;
                nu[1] = nu[1] * cos_theta + v[1] * sin_theta + e[1] * factor_2;
                nu[2] = nu[2] * cos_theta + v[2] * sin_theta + e[2] * factor_2;
                
                // Normalize for stability
                nu.Normalize();
            }
            
            // Finally, we have to compute the vertex (n-1). We need to apply only an alpha-rotation.
            
            // Compute angle alpha beteen diagonals d[i-1] and d[i] by the cosine identity.
            
            const Real cos_alpha = ( d[n-3] * d[n-3] + d[n-2] * d[n-2] - one )/( two * d[n-3] * d[n-2] );
            
            // 0 < alpha < Pi, so sin(alpha) is postive. Thus the following is safe.
            const Real sin_alpha = std::sqrt(one - cos_alpha * cos_alpha);
            
            // Cross product of nu and unit vector e.
            Cross(nu,e,v);
            
            const Real factor = Dot(nu,e) * (one-cos_alpha);
            
            // Apply Rodrigues' formula
            e[0] = e[0] * cos_alpha + v[0] * sin_alpha + nu[0] * factor;
            e[1] = e[1] * cos_alpha + v[1] * sin_alpha + nu[1] * factor;
            e[2] = e[2] * cos_alpha + v[2] * sin_alpha + nu[2] * factor;
            
            // Normalize for stability
            e.Normalize();
            
            // Compute the new vertex position.
            p[3 * (n-1) + 0] = d[n-2] * e[0];
            p[3 * (n-1) + 1] = d[n-2] * e[1];
            p[3 * (n-1) + 2] = d[n-2] * e[2];
            
            // Fill the last (redundant) vertex with 0.
            p[3 * n + 0] = Real(0);
            p[3 * n + 1] = Real(0);
            p[3 * n + 2] = Real(0);
            
            return trials;
        }
        
        /*!
         * Generates `sample_count` closed, equilateral random polygons.
         *
         * This is a port of the routine `plc_random_equilateral_closed_polygon` from the C library [plCurve](https://jasoncantarella.com/wordpress/software/plcurve) by Ted Ashton, Jason Cantarella, Harrison Chapman, and Tom Eddy.
         *
         * @param p Buffer for vertex positions; assumed to be of size `sample_count * EdgeCount() * AmbientDimension()`. Coordinates are stored in interleaved form, i.e. the coordinates of each vertex lie contiguously in memory.
         *
         * @param sample_count Number of samples to generate.
         *
         * @param thread_count Number of threads to use. Best practice is to set this to the number of performance cores on your system.
         */
        
        Int CreateRandomClosedPolygons( Real * restrict const p, const Int sample_count, const Int thread_count = 1 )
        {
            ptic(ClassName()+"::CreateRandomClosedPolygons");
            
            const Int trials = ParallelDoReduce(
                [=]( const Int thread) -> Int
                {
                    
                    const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                    const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );
                    
                    // Create a new instance of the class with its own random number generator.
                    Sampler C ( edge_count_ );
                
                    Int trials = 0;
                    
                    const Int step = AmbDim * (edge_count_ + 1);
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        trials += C.CreateRandomClosedPolygon( &p[ step * k ] );
                    }
                    
                    return trials;
                },
                AddReducer<Int,Int>(),
                Scalar::Zero<Int>,
                thread_count
            );
            
            ptoc(ClassName()+"::CreateRandomClosedPolygons");
            return trials;
        }
        
    public:
        
        /*!
         * @brief Returns the dimension of the ambient space.
         */
        
        static constexpr Int AmbientDimension()
        {
            return AmbDim;
        }
        
        /*! @brief Returns a string that identifies the class. Good for debugging and printing messages. */
        
        std::string ClassName()
        {
            return std::string("AAM::Sampler") + "<" + TypeName<Real> + "," + TypeName<Int>  + "," + random_engine.ClassName() + ","+ Tools::ToString(ProgressiveQ)+ ">";
        }
    };
    
} // namespace AAM
