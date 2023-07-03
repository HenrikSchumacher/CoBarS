#pragma once

namespace CycleSampler
{

    template<typename Real, typename Int>
    class MomentPolytopeSampler
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using SpherePoints_T = Tensor2<Real,Int>;
        using SpacePoints_T  = Tensor2<Real,Int>;
        
        static constexpr Int AmbDim = 3;
        
        using Vector_T = Tensors::Tiny::Vector<AmbDim,Real,Int>;
        
        MomentPolytopeSampler() = default;
        
        ~MomentPolytopeSampler(){}
        
        explicit MomentPolytopeSampler( const Int edge_count_ )
        :   edge_count  ( edge_count_)
        {
        }

    protected:
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
        static constexpr Real two_pi            = static_cast<Real>(2 * M_PI);
        
        const Int edge_count;
        
        std::mt19937_64 engine { std::random_device()() };
                    
        std::uniform_real_distribution<Real> dist_1 { static_cast<Real>(-1), static_cast<Real>(1) };
        
        std::uniform_real_distribution<Real> dist_2 { static_cast<Real>(0), static_cast<Real>(2 * M_PI) };
        
        static constexpr Int max_trials = 10000000;
        
    public:
        
//        const SpherePoints_T & EdgeCoordinates() const
//        {
//            return y;
//        }
//
//        void ReadEdgeCoordinates( const Real * const y_in )
//        {
//            y.Read(y_in);
//        }
//
//        void ReadEdgeCoordinates( const Real * const y_in, const Int k )
//        {
//            ReadEdgeCoordinates( &y_in[ AmbDim * edge_count * k ]);
//        }
//
//        void WriteEdgeCoordinates( Real * y_out ) const
//        {
//            y.Write(y_out);
//        }
//
//        void WriteEdgeCoordinates( Real * y_out, const Int k ) const
//        {
//            WriteEdgeCoordinates( &y_out[ AmbDim * edge_count * k ]);
//        }
        
    public:
        
        
        Int RandomClosedPolygon( mut<Real> p )
        {
            // Port of the routine "plc_random_equilateral_closed_polygon" from the C library "plCurve" by Ted Ashton, Jason Cantarella, Harrison Chapman, and Tom Eddy.
            // https://jasoncantarella.com/wordpress/software/plcurve/
            
            // We assume that user-supplied buffer p has size edge_count * AmbDim;
            
            // Storing the coordinates in interleaved form, i.e.
            // p = { x[0], y[0], z[0], x[1], y[1], z[1], x[2], y[2], z[2], ... }
            
            const Int n = edge_count;
            
            // We use the user-supplied buffer as scratch space for the diagonal lengths.
            mut<Real> d = &p[2*edge_count+1];
            
            bool rejected = true;

            Int trials = 0;
            
            // We need to find d[0], d[1], ... , d[n-2] from the moment polytope.
            // d[0], d[1], ..., d[n-2] are diagonals of a fan polyhedron if and only if the following is satisfied:
            //
            // (1)      d[0  ] = 1;
            // (2)      d[n-2] = 1;
            
            // (4)      |d[i-1] - d[i]| <= 1    for i = 1,...,n-3
            // (5)      d[i-1] + d[i] >= 1      for i = 1,...,n-3
            // (6)      d[i] >0                 for i = 1,...,n-3
            
            // (7)      |d[n-3] - d[n-2]| <= 1
            // (8)      d[n-3] + d[n-2] >= 1
            
            //Because d[n-2] = 1, the latter two reduce to
            
            // (7')     |d[n-3] - 1| <= 1
            // (8')     d[n-3] >= 0             obsolete, follows from (6) for i = n-3.
        
            //Because d[n-3] >0, the condition (4') reduces to
            
            //  (7'') d[n-3] <= 2.
            
            // This guarantees (1), (2).
            d[0  ] = one;
            d[n-2] = one;
            
            while( rejected && (trials < max_trials) )
            {
                ++trials;
                
                for( Int i = 1; i < n-2; ++i )
                {
                    // This guarantees (4):
                    d[i] = d[i-1] + dist_1( engine );
                    
                    //Check condition (5) and (6).
                    rejected = (d[i-1] + d[i] < one) || (d[i] < eps);
                    
                    if( rejected )
                    {
                        break;
                    }
                }
                
                // Check condition (7'') for last diagonal:
                rejected = rejected || ( d[n-3] > 2 );
            }
            
            if( trials > max_trials )
            {
                eprint("No success after "+ToString(max_trials)+" trials.");
                return trials;
            }
            
//            // Draw uniformly distributed angles:
//            for( Int i = 0; i < n-3; ++i )
//            {
//                theta[i] = dist_2( engine );
//            }
            
            p[0] = zero;
            p[1] = zero;
            p[2] = zero;
            p[3] = one ;
            p[4] = zero;
            p[5] = zero;
            
            // Unit vector pointing to position of currect vertex..
            Vector_T e;
            e[0] = one;
            e[1] = zero;
            e[2] = zero;
            
            // Some buffer for the cross product in Rodrigue's formula.
            Vector_T v;
            
            // The current triangle's normal.
            Vector_T nu;
            e[0] = zero;
            e[1] = zero;
            e[2] = one;
            
            for( Int i = 0; i < n-3; ++i )
            {
                // Next we compute the new unit vector that points to e by rotating e by the angle alpha about the unit normal of the triangle.
                
                // Compute angle alpha beteen diagonals d[i-1] and d[i] by the cosine identity.
                const Real cos_alpha = ( d[i] * d[i] + d[i+1] * d[i+1] - one )/( two * d[i] * d[i+1] );
                
                // 0 < alpha < Pi, so sin(alpha) is postive. Thus the following is safe.
                const Real sin_alpha = std::sqrt(one - cos_alpha*cos_alpha);
                
                // Cross product of nu and unit vector e.
                v[0] = nu[1] * e[2] - nu[2] * e[1];
                v[1] = nu[2] * e[0] - nu[0] * e[2];
                v[2] = nu[0] * e[1] - nu[1] * e[0];
                
                const Real factor = Dot(nu,e) * (one-cos_alpha);
                
                // Apply Rodrigue's formula
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
                
                // Cross product of e and nu.
                v[0] = e[1] * nu[2] - e[2] * nu[1];
                v[1] = e[2] * nu[0] - e[0] * nu[2];
                v[2] = e[0] * nu[1] - e[1] * nu[0];
                
                const Real theta_i   = dist_2( engine );
                const Real cos_theta = std::cos(theta_i);
                const Real sin_theta = std::sin(theta_i);
                const Real factor_2 = Dot(e,nu) * (one-cos_theta);
                
                // Apply Rodrigue's formula
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
            const Real sin_alpha = std::sqrt(one - cos_alpha*cos_alpha);
            
            // Cross product of nu and unit vector e.
            v[0] = nu[1] * e[2] - nu[2] * e[1];
            v[1] = nu[2] * e[0] - nu[0] * e[2];
            v[2] = nu[0] * e[1] - nu[1] * e[0];
            
            const Real factor = Dot(nu,e) * (one-cos_alpha);
            
            // Apply Rodrigue's formula
            e[0] = e[0] * cos_alpha + v[0] * sin_alpha + nu[0] * factor;
            e[1] = e[1] * cos_alpha + v[1] * sin_alpha + nu[1] * factor;
            e[2] = e[2] * cos_alpha + v[2] * sin_alpha + nu[2] * factor;
            
            // Normalize for stability
            e.Normalize();
            
            // Compute the new vertex position.
            p[3 * (n-1) + 0] = d[n-2] * e[0];
            p[3 * (n-1) + 1] = d[n-2] * e[1];
            p[3 * (n-1) + 2] = d[n-2] * e[2];
            
            return trials;
        }
        
        
        Int RandomClosedPolygons( mut<Real> p, const Int sample_count, const Int thread_count = 1 )
        {
            ptic(ClassName()+"::RandomClosedPolygons");
            
            const Int trials = ParallelDoReduce(
                [=]( const Int thread) -> Int
                {
                    
                    const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                    const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );
                    
                    // Create a new instance of the class with its own random number generator.
                    MomentPolytopeSampler C ( edge_count );
                
                    Int trials = 0;
                    
                    const Int step = AmbDim * edge_count;
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        trials += C.RandomClosedPolygon( &p[ step * k ] );
                    }
                    
                    return trials;
                },
                AddReducer<Int,Int>(),
                Scalar::Zero<Int>,
                thread_count
            );
            
            ptoc(ClassName()+"::RandomClosedPolygons");
            return trials;
        }
        
    public:
        
        static constexpr Int AmbientDimension()
        {
            return AmbDim;
        }
        
        std::string ClassName()
        {
            return std::string("MomentPolytopeSampler") + "<" + TypeName<Real> + "," + TypeName<Int> + ">";
        }
    };
    
} // namespace CycleSampler
