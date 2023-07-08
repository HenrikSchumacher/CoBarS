#pragma once

namespace CycleSampler
{
    
    template<int AmbDim, typename Real = double, typename Int = int_fast32_t, typename PRNG_T = Xoshiro256Plus,
        bool zerofy_first = true
    >
    class Sampler_vec : public SamplerBase<AmbDim,Real,Int>
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    private:
        
        using Base_T = SamplerBase<AmbDim,Real,Int>;
        
    public:
        
        using Vector_T          = typename Base_T::Vector_T;
        using SquareMatrix_T    = typename Base_T::SquareMatrix_T;
        using SymmetricMatrix_T = typename Base_T::SymmetricMatrix_T;
        
        using RandomVariable_T  = typename Base_T::RandomVariable_T;
        
        
        using Weights_T         = typename Base_T::Weights_T;
        using Setting_T         = typename Base_T::Setting_T;
        
        using SpherePoints_T    = Tiny::VectorList<AmbDim,Real,Int>;
        using SpacePoints_T     = Tiny::VectorList<AmbDim,Real,Int>;

        
        Sampler_vec() = default;
        
        ~Sampler_vec() = default;
        
        explicit Sampler_vec(
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   Base_T( settings_ )
        ,   edge_count(edge_count_)
        ,   x   ( edge_count )
        ,   y   ( edge_count )
        ,   p   ( edge_count + 1 )
        ,   r   ( edge_count, one / edge_count )
        ,   rho ( edge_count, one )
        ,   total_r_inv ( one )
        {}
        
        explicit Sampler_vec(
            ptr<Real> r_in,
            ptr<Real> rho_in,
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   Base_T( settings_ )
        ,   edge_count(edge_count_)
        ,   x   ( edge_count )
        ,   y   ( edge_count )
        ,   p   ( edge_count + 1 )
        ,   r   ( edge_count )
        ,   rho ( rho_in, edge_count )
        {
            ReadEdgeLengths(r_in);
        }
        
        
        // Copy constructor
        Sampler_vec( const Sampler_vec & other )
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
        ,   edge_space_sampling_weight(other.edge_space_sampling_weight)
        ,   edge_quotient_space_sampling_correction(other.edge_quotient_space_sampling_correction)
        ,   lambda_min(other.lambda_min)
        ,   q(other.q)
        ,   errorestimator(errorestimator)
        ,   linesearchQ(linesearchQ)
        ,   succeededQ(succeededQ)
        ,   continueQ(continueQ)
        ,   ArmijoQ(ArmijoQ)
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                random_engine[i].seed(pcg_extras::seed_seq_from<std::random_device>());
            }
        }
        
        friend void swap(Sampler_vec &A, Sampler_vec &B) noexcept
        {
            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
            using std::swap;
            swap(A.edge_count,B.edge_count);
            swap(A.settings,B.settings);
            swap(A.x,B.x);
            swap(A.y,B.y);
            swap(A.p,B.p);
            swap(A.r,B.r);
            swap(A.rho,B.rho);
            swap(A.total_r_inv,B.total_r_inv);
            swap(A.w,B.w);
            swap(A.F,B.F);
            swap(A.A,B.A);
            swap(A.u,B.u);
            swap(A.z,B.z);

            swap(A.iter,             B.iter             );
            swap(A.squared_residual, B.squared_residual );
            swap(A.residual,         B.residual         );
            swap(A.edge_space_sampling_weight,A.edge_space_sampling_weight);
            swap(A.edge_quotient_space_sampling_correction,A.edge_quotient_space_sampling_correction);
            swap(A.lambda_min,B.lambda_min);
            swap(A.q,B.q);
            swap(A.errorestimator,B.errorestimator);
            swap(A.linesearchQ,B.linesearchQ);
            swap(A.succeededQ,B.succeededQ);
            swap(A.continueQ,B.continueQ);
            swap(A.ArmijoQ,B.ArmijoQ);
        }
        
        // Copy assignment operator
        Sampler_vec & operator=(Sampler_vec other)
        {
            // copy-and-swap idiom
            // see https://stackoverflow.com/a/3279550/8248900 for details
            swap(*this, other);

            return *this;
        }

        /* Move constructor */
        Sampler_vec( Sampler_vec && other ) noexcept
        :   Sampler_vec()
        {
            swap(*this, other);
        }
        
        
    protected:
        
        const Int edge_count = 0;
        
        mutable PRNG_T random_engine [AmbDim];
        
        mutable std::normal_distribution<Real> normal_dist {zero,one};
        
        SpherePoints_T x {0};
        SpherePoints_T y {0};
        SpacePoints_T  p {0};
        
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
        
        Real edge_space_sampling_weight = 0;
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
    
    protected:
        
    public:
        
        using Base_T::Settings;
        
        virtual void Optimize() override
        {
            const Int max_iter = Settings().max_iter;
            
            iter = 0;
            
            Shift();
            
            DifferentialAndHessian_Hyperbolic();
            
            SearchDirection_Hyperbolic();
            
            // The main loop.
            while( ( iter < max_iter ) && continueQ )
            {
                ++iter;
                
                LineSearch_Hyperbolic_Potential();
                
                DifferentialAndHessian_Hyperbolic();
                
                SearchDirection_Hyperbolic();
            }
        }
        
        virtual Int EdgeCount() const override
        {
            return edge_count;
        }
        
        virtual Real EdgeSpaceSamplingWeight() const override
        {
            return edge_space_sampling_weight;
        }
        
        virtual Real EdgeQuotientSpaceSamplingCorrection() const override
        {
            return edge_quotient_space_sampling_correction;
        }
        
        virtual Real EdgeQuotientSpaceSamplingWeight() const override
        {
            return EdgeSpaceSamplingWeight() * EdgeQuotientSpaceSamplingCorrection();
        }
        
        
    protected:
        
        Real Potential()
        {
            const Real zz = Dot(z,z);
            
            const Real a = big_one + zz;
            const Real c = (big_one-zz);
            
            const Real b = one/c;
            
            Real value = 0;
            
            for( Int k = 0; k < edge_count; ++k )
            {
                Vector_T y_k ( y, k );
                
                value += r[k] * std::log( std::abs( (a - two * Dot(y_k,z) ) * b ) );
            }
            
            return value * total_r_inv;
        }
        
        
        void LineSearch_Hyperbolic_Residual()
        {
            // 2 F(0)^T.DF(0).u is the derivative of w\mapsto F(w)^T.F(w) at w = 0.
            const Real slope = two * DF.InnerProduct(F,u);
            
            Real tau = one;
            
            const Real u_norm = u.Norm();
            
            // exponential map shooting from 0 to tau * u.
            Times( tau * tanhc(tau * u_norm), u, z );

            // Shift the point z along -w to get new updated point w
            InverseShift();
            
            // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
            Shift();
            
            const Real squared_residual_at_0 = squared_residual;
            
            DifferentialAndHessian_Hyperbolic();
            
            if( linesearchQ )
            {
                Int backtrackings = 0;
                
                // Armijo condition for squared residual.
                
                ArmijoQ = squared_residual - squared_residual_at_0 - Settings().Armijo_slope_factor * tau * slope < zero;
                
                while( !ArmijoQ && (backtrackings < Settings().max_backtrackings) )
                {
                    ++backtrackings;
                    
                    // Estimate step size from quadratic fit if applicable.
                    
                    const Real tau_1 = Settings().Armijo_shrink_factor * tau;
                    const Real tau_2 = - half * Settings().Armijo_slope_factor * tau * tau * slope / ( squared_residual  - squared_residual_at_0 - tau * slope );
                    
                    tau = std::max( tau_1, tau_2 );
                    
                    Times( tau * tanhc(tau * u_norm), u, z );
                    
                    // Shift the point z along -w to get new updated point w .
                    InverseShift();
                    
                    // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
                    Shift();
                    
                    DifferentialAndHessian_Hyperbolic();
                    
                    ArmijoQ = squared_residual - squared_residual_at_0 - Settings().Armijo_slope_factor * tau * slope < 0;
                }
            }
        }
        
        void LineSearch_Hyperbolic_Potential()
        {
            Real tau = one;
            
            const Real u_norm = u.Norm();
            
            // exponential map shooting from 0 to tau * u.
            
            Times( tau * tanhc(tau * u_norm), u, z );
            
            if( linesearchQ )
            {
                //Linesearch with potential as merit function.
                
                const Real gamma = Settings().Armijo_shrink_factor;
                
                const Real sigma = Settings().Armijo_slope_factor;
                
                const Real Dphi_0 = g_factor * Dot(F,u);
                
                Int backtrackings = 0;
                
                // Compute potential and check Armijo condition.
                
                // const Real phi_0 = 0;
                
                Real phi_tau = Potential();
                
                ArmijoQ = phi_tau /*- phi_0*/ - sigma * tau * Dphi_0 < 0;
                
                
                while( !ArmijoQ && (backtrackings < Settings().max_backtrackings) )
                {
                    ++backtrackings;
                    
                    const Real tau_1 = gamma * tau;
                    
                    // Estimate step size from quadratic fit if applicable.
                    const Real tau_2 = - half * sigma * tau * tau * Dphi_0 / ( phi_tau /*- phi_0*/ - tau * Dphi_0 );
                    
                    tau = std::max( tau_1, tau_2 );
                    
                    Times( tau * tanhc(tau * u_norm), u, z );
                    
                    phi_tau = Potential();
                    
                    ArmijoQ = phi_tau  /*- phi_0*/ - sigma * tau * Dphi_0 < 0;
                }
            }
            
            // Shift the point z along -w to get new updated point w .
            InverseShift();
            
            // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
            Shift();
        }
        
        void DifferentialAndHessian_Hyperbolic()
        {
            // CAUTION: We use a different sign convention as in the paper!
            // Assemble  F = -1/2 y * r.
            // Assemble DF = nabla F + regulatization:
            // DF_{ij} = \delta_{ij} - \sum_k x_{k,i} x_{k,j} \r_k.
            
            if constexpr ( zerofy_first )
            {
                F.SetZero();
                DF.SetZero();
                
                for( Int k = 0; k < edge_count; ++k )
                {
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        const Real factor = r[k] * y[i][k];
                        
                        F[i] -= factor;
                        
                        for( Int j = i; j < AmbDim; ++j )
                        {
                            DF[i][j] -= factor * y[j][k];
                        }
                    }
                }
            }
            else
            {
                // Filling F and DF with first summand...
                {
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        const Real factor = r[0] * y[i][0];
                        
                        F[i] = - factor;
                        
                        for( Int j = i; j < AmbDim; ++j )
                        {
                            DF[i][j] = - factor * y[j][0];
                        }
                    }
                }
                
                // ... and adding-in the other summands.
                for( Int k = 1; k < edge_count; ++k )
                {
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        const Real factor = r[k] * y[i][k];
                        
                        F[i] -= factor;
                        
                        for( Int j = i; j < AmbDim; ++j )
                        {
                            DF[i][j] -= factor * y[j][k];
                        }
                    }
                }
            }
            
            // Normalize for case that the weights in r do not sum to 1.
            
            F  *= total_r_inv;
            DF *= total_r_inv;
            
            squared_residual = Dot(F,F);
            
            residual = std::sqrt( squared_residual );
            
            F *= half;
            
            // Better add the identity afterwards for precision reasons.
            for( Int i = 0; i < AmbDim; ++i )
            {
                DF[i][i] += one;
            }
        }
        
        void SearchDirection_Hyperbolic()
        {
            // Make decisions whether to continue.
            if( residual < static_cast<Real>(100.) * Settings().tolerance )
            {
                // We have to compute eigenvalue _before_ we add the regularization.
                lambda_min = DF.SmallestEigenvalue();
                q = four * residual / (lambda_min * lambda_min);
                
                if( q < one )
                {
                    //Kantorovich condition satisfied; this allows to compute an error estimator.
                    errorestimator = half * lambda_min * q;
                    //And we should deactivate line search. Otherwise, we may run into precision issues.
                    linesearchQ = false;
                    continueQ = (errorestimator > Settings().tolerance);
                    succeededQ = !continueQ;
                }
                else
                {
                    errorestimator = infty;
                    linesearchQ = Settings().Armijo_slope_factor > zero;
                    //There is no way to reduce the residual below machine epsilon. If the algorithm reaches here, the problem is probably too ill-conditioned to be solved in machine precision.
                    continueQ = residual > Settings().give_up_tolerance;
                }
            }
            else
            {
                q = big_one;
                lambda_min = eps;
                errorestimator = infty;
                linesearchQ = Settings().Armijo_slope_factor > zero;
                continueQ = residual>std::max( Settings().give_up_tolerance, Settings().tolerance );
            }
            
            const Real c = Settings().regularization * squared_residual;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                for( Int j = i; j < AmbDim; ++j )
                {
                    L[i][j] = DF[i][j] + static_cast<Real>(i==j) * c;
                }
            }
            
            L.Cholesky();
            
            L.CholeskySolve(F,u);
            
            u *= -one;
        }
        
        void Gradient_Hyperbolic()
        {
            Times(-g_factor_inv, F, u);
        }
        
        void Gradient_Planar()
        {
            // The factor 2 is here to reproduce the Abikoff-Ye algorithm (in the absence of linesearch.)
            Times(-two, F, u);
        }
        
        void InverseShift()
        {
            const Real ww  = Dot(w,w);
            const Real wz2 = Dot(w,z) * two;
            const Real zz  = Dot(z,z);
            
            const Real a = one - ww;
            const Real b = one + zz + wz2;
            const Real c = big_one + wz2 + ww * zz;
            const Real d = one / c;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                w[i] = ( a * z[i] + b * w[i] ) * d;
            }
        }
        
        void Shift()
        {
            // Shifts all entries of x along w and writes the results to y.
            
            const Real ww = Dot(w,w);
            const Real one_minus_ww = big_one - ww;
            const Real one_plus_ww  = big_one + ww;
            
            if( ww <= norm_threshold )
            {
                for( Int k = 0; k < edge_count; ++k )
                {
                    Vector_T x_k ( x, k );
                                         
                    const Real wx2 = two * Dot(w,x_k);
                    
                    const Real denom = one / ( one_plus_ww - wx2 );
                    
                    const Real wx2_minus_2 = wx2 - two;
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        y[i][k] = (one_minus_ww * x_k[i] + wx2_minus_2 * w[i]) * denom;
                    }
                }
            }
            else
            {
                // If w lies close to the boundary of the ball, then normalizing the output is a good idea.
                
                for( Int k = 0; k < edge_count; ++k )
                {
                    Vector_T x_k ( x, k );
                    
                    const Real wx2 = two * Dot(w,x_k);
                    
                    const Real denom = one / ( one_plus_ww - wx2 );
                    
                    const Real wx2_minus_2 = wx2 - two;
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        x_k[i] = (one_minus_ww * x_k[i] + wx2_minus_2 * w[i]) * denom;
                    }
                    
                    x_k.Normalize();
                    
                    x_k.Write( y, k );
                }
            }
        }
        
        void ComputeEdgeSpaceSamplingWeight()
        {
            // Shifts all entries of x along y and writes the results to y.
            // Mind that x and y are stored in SoA fashion, i.e., as matrix of size AmbDim x point_count.
            
            SquareMatrix_T cbar  (zero);
            SquareMatrix_T gamma (zero);
            
            Real prod = one;
            
            const Real ww = Dot(w,w);
            
            const Real one_plus_ww = big_one + ww;
            
            for( Int k = 0; k < edge_count; ++k )
            {
                Vector_T y_k ( y, k );
                
                const Real wy = Dot(w,y_k);
                
                const Real factor = one_plus_ww + two * wy;
                
                // Multiplying by one_plus_ww_inv so that prod does not grow so quickly.
                prod *= factor;
                
                const Real r_k = r[k];
                const Real r_over_rho_k = r_k / rho[k];
                const Real r_over_rho_k_squared = r_over_rho_k * r_over_rho_k;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    for( Int j = 0; j < AmbDim; ++j )
                    {
                        const Real scratch = (static_cast<Real>(i==j) - y_k[i] *  y_k[j] );
                        
                        gamma[i][j] += r_over_rho_k_squared * scratch;
                        
                        cbar [i][j] += r_k * scratch;
                    }
                }
            }
            
            // We can simply absorb the factor std::pow(2/(one_minus_ww),d) into the function chi.
            //  cbar *= static_cast<Real>(2)/(one_minus_ww);
            
            edge_space_sampling_weight = Power(prod, static_cast<Int>(AmbDim-1)) * sqrt(gamma.Det()) / cbar.Det();
        }
        
        void ComputeEdgeQuotientSpaceSamplingCorrection()
        {
            if constexpr ( AmbDim == 2)
            {
                edge_quotient_space_sampling_correction = one;
                return;
            }
            
            Tiny::SelfAdjointMatrix<AmbDim, Real, Int> Sigma;
            
            // We fill only the upper triangle of Sigma, because that's the only thing that the function Eigenvalues needs.
          
            
            if constexpr ( zerofy_first )
            {
                Sigma.SetZero();
                
                for( Int k = 0; k < edge_count; ++k )
                {
                    const Real rho_squared = rho[k] * rho[k];
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        const Real factor = rho_squared * y[i][k];
                        
                        for( Int j = i; j < AmbDim; ++j )
                        {
                            Sigma[i][j] += factor * y[j][k];
                        }
                    }
                }
            }
            else
            {
                // Overwrite for k = 0.
                {
                    const Real rho_squared = rho[0] * rho[0];
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        const Real factor = rho_squared * y[i][0];

                        for( Int j = i; j < AmbDim; ++j )
                        {
                            Sigma[i][j] = factor * y[j][0];
                        }
                    }
                }

                // Now we add-in the other entries.
                for( Int k = 1; k < edge_count; ++k )
                {
                    const Real rho_squared = rho[k] * rho[k];
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        const Real factor = rho_squared * y[i][k];

                        for( Int j = i; j < AmbDim; ++j )
                        {
                            Sigma[i][j] += factor * y[j][k];
                        }
                    }
                }
            }
            
            
            
            if constexpr ( AmbDim == 3)
            {
                // Exploiting that
                //      (lambda[0] + lambda[1]) * (lambda[0] + lambda[2]) * (lambda[1] + lambda[2])
                //      =
                //      ( tr(Sigma*Sigma) - tr(Sigma)*tr(Sigma) ) *  tr(Sigma)/2 - det(Sigma)
                //  Thus, it can be expressed by as third-order polynomial in the entries of the matrix.
                
                const Real S_00 = Sigma[0][0] * Sigma[0][0];
                const Real S_11 = Sigma[1][1] * Sigma[1][1];
                const Real S_22 = Sigma[2][2] * Sigma[2][2];
                
                const Real S_10 = Sigma[0][1]*Sigma[0][1];
                const Real S_20 = Sigma[0][2]*Sigma[0][2];
                const Real S_21 = Sigma[1][2]*Sigma[1][2];
                
                const Real det = std::abs(
                      Sigma[0][0] * ( S_11 + S_22 - S_10 - S_20 )
                    + Sigma[1][1] * ( S_00 + S_22 - S_10 - S_21 )
                    + Sigma[2][2] * ( S_00 + S_11 - S_20 - S_21 )
                    + two * (Sigma[0][0]*Sigma[1][1]*Sigma[2][2] - Sigma[0][1]*Sigma[0][2]*Sigma[1][2])
                );
                edge_quotient_space_sampling_correction = one / std::sqrt(det);
            }
            else
            {
                Tiny::Vector<AmbDim,Real,Int> lambda;
                
                // Compute eigenvalues by QR algorithm.
                
                Sigma.Eigenvalues(lambda);
                
                Real det = one;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    for( Int j = i+1; j < AmbDim; ++j )
                    {
                        det *= (lambda(i)+lambda(j));
                    }
                }
                
                edge_quotient_space_sampling_correction = one / std::sqrt(det);
            }
        }
        
    public:
        
//        virtual const SpherePoints_T & InitialEdgeCoordinates() const override
//        {
//            return x;
//        }
        
        virtual Real InitialEdgeCoordinates( const Int k, const Int i ) const override
        {
            return x[i][k];
        }
        
        virtual Vector_T InitialEdgeCoordinates( const Int k ) const override
        {
            return Vector_T ( x, k );
        }
        
        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, bool normalize = true ) override
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
        
        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, const Int k, bool normalize = true ) override
        {
            ReadInitialEdgeCoordinates( &x_in[ AmbDim * edge_count * k], normalize );
        }
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out ) const override
        {
            x.Write( x_out );
        }
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out, const Int k ) const override
        {
            WriteInitialEdgeCoordinates( &x_out[ AmbDim * edge_count * k ]);
        }

        virtual void RandomizeInitialEdgeCoordinates() override
        {
            for( Int k = 0; k < edge_count; ++k )
            {
                Vector_T x_k;

                for( Int i = 0; i < AmbDim; ++i )
                {
                    x_k[i] = normal_dist( random_engine[i] );
                }

                x_k.Normalize();

                x_k.Write( x, k );
            }
        }
        
//        virtual const SpherePoints_T & EdgeCoordinates() const override
//        {
//            return y;
//        }
        
        virtual Real EdgeCoordinates( const Int k, const Int i ) const override
        {
            return y[i][k];
        }
        
        virtual Vector_T EdgeCoordinates( const Int k ) const override
        {
            return Vector_T ( y, k );
        }
        
        virtual void ReadEdgeCoordinates( ptr<Real> y_in ) override
        {
            y.Read( y_in );
        }
        
        virtual void ReadEdgeCoordinates( ptr<Real> y_in, const Int k ) override
        {
            ReadEdgeCoordinates( &y_in[ AmbDim * edge_count * k ]);
        }
        
        virtual void WriteEdgeCoordinates( mut<Real> y_out ) const override
        {
            y.Write( y_out );
        }
        
        virtual void WriteEdgeCoordinates( Real * y_out, const Int k ) const override
        {
            WriteEdgeCoordinates( &y_out[ AmbDim * edge_count * k ]);
        }
        
        
        
//        virtual const SpacePoints_T & SpaceCoordinates() const override
//        {
//            return p;
//        }
//        
        virtual Real SpaceCoordinates( const Int k, const Int i ) const override
        {
            return p[i][k];
        }
        
        virtual Vector_T SpaceCoordinates( const Int k ) const override
        {
            return Vector_T ( p, k );
        }
        
        virtual void WriteSpaceCoordinates( Real * p_out ) const override
        {
            p.Write( p_out );
        }
        
        virtual void WriteSpaceCoordinates( Real * p_out, const Int k ) const override
        {
            WriteSpaceCoordinates( &p_out[ (edge_count+1) * AmbDim * k ]);
        }
        
        virtual void ComputeSpaceCoordinates() override
        {
            //Caution: This gives only have the weight to the end vertices of the chain.
            //Thus this is only really the barycenter, if the chain is closed!
            
            Vector_T barycenter        (zero);
            Vector_T point_accumulator (zero);
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Real r_k = r[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real offset = r_k * y[i][k];
                    
                    barycenter[i] += (point_accumulator[i] + half * offset);
                    
                    point_accumulator[i] += offset;
                }
            }
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                p[i][0] = -barycenter[i] / edge_count;
            }
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Real r_k = r[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    p[i][k+1] = p[i][k] + r_k * y[i][k];
                }
            }
        }
        
        
        virtual const Weights_T & EdgeLengths() const override
        {
            return r;
        }
        
        virtual void ReadEdgeLengths( const Real * const r_in ) override
        {
            r.Read(r_in);
            
            total_r_inv = one / r.Total();
        }
        
        
        virtual const Weights_T & Rho() const override
        {
            return rho;
        }
        
        virtual void ReadRho( const Real * const rho_in ) override
        {
            rho.Read(rho_in);
        }
        
        
        virtual void ComputeShiftVector() override
        {
            w.SetZero();
            
            if constexpr ( zerofy_first )
            {
                for( Int k = 0; k < edge_count; ++k )
                {
                    const Real r_k = r[k];
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        w[i] += x[i][k] * r_k;
                    }
                }
            }
            else
            {
                // Overwrite by first summand.
                {
                    const Real r_k = r[0];

                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        w[i] = x[i][0] * r_k;
                    }
                }

                // Add-in the others.
                for( Int k = 1; k < edge_count; ++k )
                {
                    const Real r_k = r[k];

                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        w[i] += x[i][k] * r_k;
                    }
                }
            }

            
            // Normalize in that case that r does not sum up to 1.
            w *= total_r_inv;
        }
        
        virtual void ReadShiftVector( const Real * const w_in ) override
        {
            w.Read(w_in);
            
            // Use Euclidean barycenter as initial guess if the supplied initial guess does not make sense.
            if( Dot(w,w) > small_one )
            {
                ComputeShiftVector();
            }
        }
        
        virtual void ReadShiftVector( const Real * const w_in, const Int k ) override
        {
            ReadShiftVector( &w_in[ AmbDim * k ] );
        }
        
        virtual void WriteShiftVector( Real * w_out ) const override
        {
            w.Write( w_out );
        }
        
        virtual void WriteShiftVector( Real * w_out, const Int k ) const override
        {
            w.Write( &w_out[ AmbDim * k] );
        }
        
        virtual const Vector_T & ShiftVector() const override
        {
            return w;
        }
        
        
        virtual Real Residual() const override
        {
            return residual;
        }
        
        virtual Real ErrorEstimator() const override
        {
            return errorestimator;
        }
        
        virtual Int IterationCount() const override
        {
            return iter;
        }
        
        virtual Int MaxIterationCount() const override
        {
            return Settings().max_iter;
        }
        
        
    public:
        
        virtual void OptimizeBatch(
            ptr<Real>  x_in,
            mut<Real>  w_out,
            mut<Real>  y_out,
            const Int  sample_count,
            const Int  thread_count = 1,
            const bool normalize = true
        ) override
        {
            ptic(ClassName()+"::OptimizeBatch");
            
            ParallelDo(
                [&,this]( const Int thread )
                {
                    const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                    const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );

                    Sampler_vec S( edge_count, Settings() );

                    S.ReadEdgeLengths( EdgeLengths().data() );

                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        S.ReadInitialEdgeCoordinates( x_in, k, normalize );

                        S.ComputeShiftVector();

                        S.Optimize();

                        S.WriteShiftVector( w_out, k );

                        S.WriteEdgeCoordinates( y_out, k );
                    }
                },
                thread_count
            );
            
            ptoc(ClassName()+"::OptimizeBatch");
        }
        
        virtual void RandomClosedPolygons(
            mut<Real> x_out,
            mut<Real> w_out,
            mut<Real> y_out,
            mut<Real> K_edge_space,
            mut<Real> K_edge_quotient_space,
            const Int sample_count,
            const Int thread_count = 1
        ) const override
        {
            ptic(ClassName()+"::RandomClosedPolygons");
            
            ParallelDo(
                [&,this]( const Int thread )
                {
                    Time start = Clock::now();
                    
                    const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                    const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );

                    Sampler_vec S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );

                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        S.RandomizeInitialEdgeCoordinates();
                        
                        S.WriteInitialEdgeCoordinates(x_out, k);
                        
                        S.ComputeShiftVector();
                        
                        S.Optimize();
                        
                        S.WriteShiftVector(w_out,k);
                        
                        S.WriteEdgeCoordinates(y_out,k);
                        
                        S.ComputeEdgeSpaceSamplingWeight();
                        
                        S.ComputeEdgeQuotientSpaceSamplingCorrection();
                        
                        K_edge_space[k] = S.EdgeSpaceSamplingWeight();
                        
                        K_edge_quotient_space[k] = S.EdgeQuotientSpaceSamplingWeight();
                    }
                    
                    Time stop = Clock::now();
                    
                    logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Duration(start, stop) ) + "." );
                    
                },
                thread_count
            );
            
            ptoc(ClassName()+"::RandomClosedPolygons");
        }
        
        virtual void Sample(
            mut<Real> sampled_values,
            mut<Real> edge_space_sampling_weights,
            mut<Real> edge_quotient_space_sampling_weights,
            std::shared_ptr<RandomVariable_T> & F_,
            const Int sample_count,
            const Int thread_count = 1
        ) const override
        {
            // This function creates sample for the random variable F and records the sampling weights, so that this weighted data can be processed elsewhere.
            
            // The generated polygons are discarded immediately after evaluating the random variables on them.
            
            // sampled_values is expected to be an array of size at least sample_count;
            // edge_space_sampling_weights is expected to be an array of size at least sample_count -- or a nullptr.
            // edge_quotient_space_sampling_weights is expected to be an array of size at least sample_count  -- or a nullptr.
            
            ptic(ClassName()+"::Sample");
            

            if( edge_space_sampling_weights != nullptr )
            {
                if( edge_quotient_space_sampling_weights != nullptr )
                {
                    sample_one<true,true>(
                        sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                        F_, sample_count, thread_count
                    );
                }
                else
                {
                    sample_one<true,false>(
                        sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                        F_, sample_count, thread_count
                    );
                }
            }
            else
            {
                if( edge_quotient_space_sampling_weights != nullptr )
                {
                    sample_one<false,true>(
                        sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                        F_, sample_count, thread_count
                    );
                }
                else
                {
                    sample_one<false,false>(
                        sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                        F_, sample_count, thread_count
                    );
                }
            }
            
            ptoc(ClassName()+"::Sample");
        }
        
    private:
        
        template<bool edge_space_flag, bool quotient_space_flag>
        void sample_one(
            mut<Real> sampled_values,
            mut<Real> edge_space_sampling_weights,
            mut<Real> edge_quotient_space_sampling_weights,
            std::shared_ptr<RandomVariable_T> & F_,
            const Int sample_count,
            const Int thread_count = 1
        ) const
        {
            ParallelDo(
                [&,this]( const Int thread )
                {
                    const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                    const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );

                    // For every thread create a copy of the current Sampler_vec object.
                    Sampler_vec S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                    
                    // Make a copy the random variable (it might have some state!).
                    std::shared_ptr<RandomVariable_T> F_ptr = F_->Clone();
                    
                    RandomVariable_T & F_loc = *F_ptr;
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        S.RandomizeInitialEdgeCoordinates();

                        S.ComputeShiftVector();

                        S.Optimize();

                        S.ComputeSpaceCoordinates();

                        if constexpr ( edge_space_flag || quotient_space_flag )
                        {
                            S.ComputeEdgeSpaceSamplingWeight();
                        }
                            
                        if constexpr ( quotient_space_flag )
                        {
                            S.ComputeEdgeQuotientSpaceSamplingCorrection();
                        }
                        
                        if constexpr ( edge_space_flag )
                        {
                            edge_space_sampling_weights[k] = S.EdgeSpaceSamplingWeight();
                        }
                        
                        if constexpr ( quotient_space_flag )
                        {
                            edge_quotient_space_sampling_weights[k] =  S.EdgeQuotientSpaceSamplingWeight();
                        }

                        sampled_values[k] = F_loc(S);
                    }
                },
                thread_count
            );
        }
        
    public:
        
        virtual void Sample(
            mut<Real> sampled_values,
            mut<Real> edge_space_sampling_weights,
            mut<Real> edge_quotient_space_sampling_weights,
            const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
            const Int sample_count,
            const Int  thread_count = 1
        ) const override
        {
            // This function creates sample for the random variables in the list F_list_ and records the sampling weights, so that this weighted data can be processed elsewhere.
            
            // The generated polygons are discarded immediately after evaluating the random variables on them.
            
            // sampled_values is expected to be an array of size at least sample_count * fun_count;
            // edge_space_sampling_weights is expected to be an array of size at least sample_count.
            // edge_quotient_space_sampling_weights is expected to be an array of size at least sample_count.
            
            ptic(ClassName()+"::Sample (batch)");

            if( edge_space_sampling_weights != nullptr )
            {
                if( edge_quotient_space_sampling_weights != nullptr )
                {
                    sample_many<true,true>(
                        sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                        F_list_, sample_count, thread_count
                    );
                }
                else
                {
                    sample_many<true,false>(
                        sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                        F_list_, sample_count, thread_count
                    );
                }
            }
            else
            {
                if( edge_quotient_space_sampling_weights != nullptr )
                {
                    sample_many<false,true>(
                        sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                        F_list_, sample_count, thread_count
                    );
                }
                else
                {
                    sample_many<false,false>(
                        sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                        F_list_, sample_count, thread_count
                    );
                }
            }
            
            ptoc(ClassName()+"::Sample (batch)");
        }
        

        template<bool edge_space_flag, bool quotient_space_flag>
        void sample_many(
            mut<Real> sampled_values,
            mut<Real> edge_space_sampling_weights,
            mut<Real> edge_quotient_space_sampling_weights,
            const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
            const Int sample_count,
            const Int  thread_count = 1
        ) const
        {
            const Int fun_count = static_cast<Int>(F_list_.size());

            ParallelDo(
                [&,this]( const Int thread )
                {
                    const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                    const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );

                    // For every thread create a copy of the current Sampler_vec object.
                    Sampler_vec S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                    
                    // Make also copys of all the random variables (they might have some state!).
                    std::vector< std::shared_ptr<RandomVariable_T> > F_list;
                    for( Int i = 0; i < fun_count; ++ i )
                    {
                        F_list.push_back(
                            std::shared_ptr<RandomVariable_T>( F_list_[i]->Clone() )
                        );
                    }
                    
                    for( Int k = k_begin; k < k_end; ++k )
                    {
                        S.RandomizeInitialEdgeCoordinates();
                        
                        S.ComputeShiftVector();
                        
                        S.Optimize();
                        
                        S.ComputeSpaceCoordinates();
                        
                        if constexpr ( edge_space_flag || quotient_space_flag )
                        {
                            S.ComputeEdgeSpaceSamplingWeight();
                        }
                        
                        if constexpr ( quotient_space_flag )
                        {
                            S.ComputeEdgeQuotientSpaceSamplingCorrection();
                        }
                        
                        if constexpr ( edge_space_flag )
                        {
                            edge_space_sampling_weights[k] = S.EdgeSpaceSamplingWeight();
                        }
                        
                        if constexpr ( quotient_space_flag )
                        {
                            edge_quotient_space_sampling_weights[k] =  S.EdgeQuotientSpaceSamplingWeight();
                        }
                        
                        for( Int i = 0; i < fun_count; ++i )
                        {
                            sampled_values[k * fun_count + i] = (*F_list[i])(S);
                        }
                    }
                },
                thread_count
            );
        }
        

        virtual void SampleCompressed(
            mut<Real> bins_out,
            const Int bin_count_,
            mut<Real> moments_out,
            const Int moment_count_,
            ptr<Real> ranges,
            const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
            const Int sample_count,
            const Int thread_count = 1
        ) const override
        {
            // This function does the sampling, but computes moments and binning on the fly, so that the sampled data can be discarded immediately.
            
            // moments: A 3D-array of size 3 x fun_count x bin_count. Entry moments(i,j,k) will store the sampled weighted k-th moment of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i (see above).
            // ranges: Specify the range for binning: For j-th function in F_list, the range from ranges(j,0) to ranges(j,1) will be devided into bin_count bins. The user is supposed to provide meaningful ranges. Some rough guess might be obtained by calling the random variables on the prepared Sampler_T C.
            
            ptic(ClassName()+"SampleCompressed");
            
            const Int fun_count = static_cast<Int>(F_list_.size());
            
            const Int moment_count = std::max( static_cast<Int>(3), moment_count_ );
            
            const Int bin_count = std::max( bin_count_, static_cast<Int>(1) );
            
            valprint( "dimension   ", AmbDim       );
            valprint( "edge_count  ", edge_count   );
            valprint( "sample_count", sample_count );
            valprint( "fun_count   ", fun_count    );
            valprint( "bin_count   ", bin_count    );
            valprint( "moment_count", moment_count );
            valprint( "thread_count", thread_count );
            
            
            Tensor3<Real,Int> bins_global   ( bins_out,    3, fun_count, bin_count    );
            Tensor3<Real,Int> moments_global( moments_out, 3, fun_count, moment_count );
            Tensor1<Real,Int> factor        (                 fun_count               );
            
            print("Sampling (compressed) the following random variables:");
            for( Int i = 0; i < fun_count; ++ i )
            {
                const size_t i_ = static_cast<size_t>(i);
                factor(i) = static_cast<Real>(bin_count) / ( ranges[2*i+1] - ranges[2*i+0] );

                print("    " + F_list_[i_]->Tag());
            }

            const Int lower = static_cast<Int>(0);
            const Int upper = static_cast<Int>(bin_count-1);
            
            
            std::mutex mutex;
            
            ParallelDo(
                [&,this]( const Int thread )
                {
                    Time start = Clock::now();
                    
                    const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                    const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );
                    
                    const Int repetitions = k_end - k_begin;
                    
                    Sampler_vec S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                    
                    std::vector< std::shared_ptr<RandomVariable_T> > F_list;
                    
                    for( Int i = 0; i < fun_count; ++ i )
                    {
                        F_list.push_back( F_list_[i]->Clone() );
                    }
                    
                    Tensor3<Real,Int> bins_local   ( 3, fun_count, bin_count,    zero );
                    Tensor3<Real,Int> moments_local( 3, fun_count, moment_count, zero );
                    
                    for( Int k = 0; k < repetitions; ++k )
                    {
                        S.RandomizeInitialEdgeCoordinates();

                        S.ComputeShiftVector();

                        S.Optimize();

                        S.ComputeSpaceCoordinates();
                        
                        S.ComputeEdgeSpaceSamplingWeight();
                        
                        S.ComputeEdgeQuotientSpaceSamplingCorrection();
                        
                        const Real K = S.EdgeSpaceSamplingWeight();

                        const Real K_quot = S.EdgeQuotientSpaceSamplingWeight();

                        for( Int i = 0; i < 1; ++i )
                        {
                            const Real val = (*F_list[i])(S);

                            Real values [3] = { one, K, K_quot };

                            const Int bin_idx = static_cast<Int>(
                                std::floor( factor[i] * (val - ranges[2*i]) )
                            );

                            if( (bin_idx <= upper) && (bin_idx >= lower) )
                            {
                                bins_local(0,i,bin_idx) += one;
                                bins_local(1,i,bin_idx) += K;
                                bins_local(2,i,bin_idx) += K_quot;
                            }

                            moments_local(0,i,0) += values[0];
                            moments_local(1,i,0) += values[1];
                            moments_local(2,i,0) += values[2];

                            for( Int j = 1; j < moment_count; ++j )
                            {
                                values[0] *= val;
                                values[1] *= val;
                                values[2] *= val;
                                moments_local(0,i,j) += values[0];
                                moments_local(1,i,j) += values[1];
                                moments_local(2,i,j) += values[2];
                            }
                        }
                    }
                    
                    {
                        const std::lock_guard<std::mutex> lock ( mutex );
                        
                        add_to_buffer<VarSize,Sequential>(
                            bins_local.data(), bins_global.data(), 3 * fun_count * bin_count
                        );
                        
                        add_to_buffer<VarSize,Sequential>(
                            moments_local.data(), moments_global.data(), 3 * fun_count * moment_count
                        );
                    }
                    
                    Time stop = Clock::now();
                 
                    logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Duration(start, stop) ) + "." );
                    
                },
                thread_count
            );
            
            bins_global.Write( bins_out );
            moments_global.Write( moments_out );
            
            ptoc(ClassName()+"::SampleCompressed");
        }

        
        
    private:
        
        Real tanhc( const Real t ) const
        {
            // Computes tanh(t)/t in a stable way by using a Padé approximation around t = 0.
            constexpr Real a0 = Scalar::One<Real>;
            constexpr Real a1 = Scalar::Frac<Real>(7,51);
            constexpr Real a2 = Scalar::Frac<Real>(1,255);
            constexpr Real a3 = Scalar::Frac<Real>(2,69615);
            constexpr Real a4 = Scalar::Frac<Real>(1,34459425);
            
            constexpr Real b0 = Scalar::One<Real>;
            constexpr Real b1 = Scalar::Frac<Real>(8,17);
            constexpr Real b2 = Scalar::Frac<Real>(7,255);
            constexpr Real b3 = Scalar::Frac<Real>(4,9945);
            constexpr Real b4 = Scalar::Frac<Real>(1,765765);
            
            const Real t2 = t * t;
            
            const Real result = ( t2 <= one )
            ? (
                a0 + t2 * (a1 + t2 * (a2 + t2 * (a3 + t2 * a4)))
            )/(
                b0 + t2 * (b1 + t2 * (b2 + t2 * (b3 + t2 * b4)))
            )
            : ( t2 <= static_cast<Real>(7) ) ? std::tanh(t)/t : one/std::abs(t);
            
            return result;
        }
        
    public:
        
        
        std::string ClassName() const override
        {
            return std::string("Sampler_vec") + "<" + ToString(AmbDim) + "," + TypeName<Real> + "," + TypeName<Int>  + "," + random_engine[0].ClassName() + ">";
        }
        
    }; // class Sampler_vec
    
} // namespace CycleSampler
