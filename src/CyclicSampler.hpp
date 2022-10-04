#pragma  once

namespace CyclicSampler {

#define CLASS CyclicSampler
#define BASE  CyclicSamplerBase<Real,Int>
    
    template<int AmbDim, typename Real = double, typename Int = long long>
    class CLASS : public BASE
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using Vector_T          = SmallVector<AmbDim,Real,Int>;
        using SymmetricMatrix_T = SmallSymmetricMatrix<AmbDim,Real,Int>;
        
        using HyperbolicPoint_T = Vector_T;
        
        using SpherePoints_T = typename BASE::SpherePoints_T;
        using SpacePoints_T  = typename BASE::SpacePoints_T;
        using Weights_T      = typename BASE::Weights_T;
        using Setting_T      = typename BASE::Setting_T;
        
        
        using BASE::settings;
        using BASE::edge_count;
        using BASE::Settings;
        using BASE::random_engine;
        using BASE::normal_dist;
        
        CLASS() : BASE() {}
        
        virtual ~CLASS(){}
        
        explicit CLASS(
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   BASE( edge_count_, settings_ )
        {
            x = SpherePoints_T( edge_count,     AmbDim );
            y = SpherePoints_T( edge_count,     AmbDim );
            p = SpherePoints_T( edge_count + 1, AmbDim );
            
            omega = Weights_T ( edge_count, one / edge_count );
            rho   = Weights_T ( edge_count, one );
            
            total_omega_inv = one;
        }
        
        explicit CLASS(
            const Real * restrict const omega_in,
            const Real * restrict const rho_in,
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   BASE( edge_count_, settings_ )
        {
            x = SpherePoints_T( edge_count,     AmbDim );
            y = SpherePoints_T( edge_count,     AmbDim );
            p = SpherePoints_T( edge_count + 1, AmbDim );
            
            omega = Weights_T( edge_count );
            rho   = Weights_T( edge_count );
            
            ReadOmega(omega_in);
            ReadRho(rho_in);
        }
        

        
//        // Copy constructor
//        CLASS( const CLASS & other )
//        :   BASE( other.edge_count, other.settings )
//        ,   x(other.x)
//        ,   y(other.y)
//        ,   p(other.p)
//        ,   omega(other.omega)
//        ,   rho(other.rho)
//        ,   total_omega_inv(other.total_omega_inv)
//        ,   w(other.w)
//        ,   F(other.F)
//        ,   DF(other.DF)
//        ,   L(other.L)
//        ,   u(other.u)
//        ,   iter(other.iter)
//        ,   residual(other.residual)
//        {}

//        
//        friend void swap(CLASS &A, CLASS &B) noexcept
//        {
//            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
//            using std::swap;
//
//            swap(A.edge_count,B.edge_count);
//            swap(A.x,B.x);
//            swap(A.y,B.y);
//            swap(A.p,B.p);
//            swap(A.omega,B.omega);
//            swap(A.rho,B.rho);
//            swap(A.total_omega_inv,B.total_omega_inv);
//            
//            std::swap_ranges(&A.w[0],    &A.w[AmbDim],            &B.w[0]   );
//            std::swap_ranges(&A.z[0],    &A.z[AmbDim],            &B.z[0]   );
//            std::swap_ranges(&A.F[0],    &A.F[AmbDim],            &B.F[0]   );
//            std::swap_ranges(&A.A(0,0),  &A.A(AmbDim-1,AmbDim),   &B.A(0,0) );
//            std::swap_ranges(&A.u[0],    &A.u[AmbDim],            &B.u[0]   );
//            
//            swap(A.settings,       B.settings     );
//            swap(A.iter,           B.iter          );
//            swap(A.residual,       B.residual      );
//        }
//        
//        // Copy assignment operator
//        CLASS & operator=(CLASS other)
//        {
//            // copy-and-swap idiom
//            // see https://stackoverflow.com/a/3279550/8248900 for details
//            swap(*this, other);
//            
//            return *this;
//        }
//        
//        /* Move constructor */
//        CLASS( CLASS && other ) noexcept
//        :   CLASS()
//        {
//            swap(*this, other);
//        }                                                                           
        

    protected:
        
        SpherePoints_T x {0,AmbDim};
        SpherePoints_T y {0,AmbDim};
        
        mutable SpacePoints_T p {0,AmbDim};
        
        Weights_T omega {0};
        Weights_T rho   {0};
        
        Real total_omega_inv = one;
        
        Vector_T w;           // current point in hyperbolic space.
        Vector_T F;           // right hand side of Newton iteration.
        SymmetricMatrix_T DF; // nabla F(0) with respect to measure ys
        SymmetricMatrix_T L;  // storing Cholesky factor.
        Vector_T u;           // update direction

        Vector_T z;           // Multiple purpose buffer.
        
        ShiftMap<AmbDim,Real,Int> S;
        
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
        static constexpr Real norm_threshold    = 0.99 * 0.99 + 16 * eps;
        static constexpr Real two_pi            = static_cast<Real>(2 * M_PI);
    
    public:
        
        virtual void Optimize() override
        {
            const Int max_iter = settings.max_iter;
            
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
        
        
        void ComputeEdgeSpaceSamplingWeight()
        {
            edge_space_sampling_weight =  S.EdgeSpaceSamplingWeight(
                x.data(), w.data(), y.data(), omega.data(), rho.data(), edge_count
            );
        }
        
        virtual Real EdgeSpaceSamplingWeight() const override
        {
            return edge_space_sampling_weight;
        }
        
        void ComputeEdgeQuotientSpaceSamplingCorrection()
        {
            if( AmbDim == 2)
            {
                edge_quotient_space_sampling_correction = one;
                return;
            }
            
            Eigen::Matrix<Real,AmbDim,AmbDim> Sigma;
            
            // We fill only the lower triangle of Sigma, because that's the only thing that Eigen' selfadjoint eigensolver needs.
            // Recall that Eigen matrices are column-major by default.
            
            {
                const Real rho_squared = rho[0] * rho[0];
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = rho_squared * y(0,i);
                    
                    for( Int j = i; j < AmbDim; ++j )
                    {
                        Sigma(j,i) = factor * y(0,j);
                    }
                }
            }
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Real rho_squared = rho[k] * rho[k];
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = rho_squared * y(k,i);
                    
                    for( Int j = i; j < AmbDim; ++j )
                    {
                        Sigma(j,i) += factor * y(k,j);
                    }
                }
            }
            
            // Eigen needs only the lower triangular part. So need not symmetrize.
            
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                for( Int j = 0; j < i; ++j )
//                {
//                    Sigma(j,i) = Sigma(i,j);
//                }
//            }

            if( AmbDim == 3)
            {
                // Exploiting that
                //      (lambda[0] + lambda[1]) * (lambda[0] + lambda[2]) * (lambda[1] + lambda[2])
                //      =
                //      ( tr(Sigma*Sigma) - tr(Sigma)*tr(Sigma) ) *  tr(Sigma)/2 - det(Sigma)
                //  Thus, it can be expressed by as third-order polynomial in the entries of the matrix.
                
                const Real S_00 = Sigma(0,0)*Sigma(0,0);
                const Real S_11 = Sigma(1,1)*Sigma(1,1);
                const Real S_22 = Sigma(2,2)*Sigma(2,2);
                
                const Real S_10 = Sigma(1,0)*Sigma(1,0);
                const Real S_20 = Sigma(2,0)*Sigma(2,0);
                const Real S_21 = Sigma(2,1)*Sigma(2,1);
                
                const Real det = std::abs(
                      Sigma(0,0) * ( S_11 + S_22 - S_10 - S_20 )
                    + Sigma(1,1) * ( S_00 + S_22 - S_10 - S_21 )
                    + Sigma(2,2) * ( S_00 + S_11 - S_20 - S_21 )
                    + two * (Sigma(0,0)*Sigma(1,1)*Sigma(2,2) - Sigma(1,0)*Sigma(2,0)*Sigma(2,1))
                );
                edge_quotient_space_sampling_correction = one / std::sqrt(det);
                return;
            }
            
            
            Eigen::SelfAdjointEigenSolver< Eigen::Matrix<Real,AmbDim,AmbDim> > eigs;
            
            eigs.compute(Sigma);

            auto & lambda = eigs.eigenvalues();

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

        virtual Real EdgeQuotientSpaceSamplingCorrection() const override
        {
            return edge_quotient_space_sampling_correction;
        }
        
        virtual Real EdgeQuotientSpaceSamplingWeight() const override
        {
            return EdgeSpaceSamplingWeight() * EdgeQuotientSpaceSamplingCorrection();
        }
        
        
    protected:
        
        Real Potential( const Vector_T & z )
        {
            Real value = 0;
            
            const Real zz = Dot(z,z);
            
            const Real a = big_one + zz;
            const Real c = (big_one-zz);
            
            const Real b = one/c;
            
            for( Int k = 0; k < edge_count; ++k )
            {
                Real yz2 = y(k,0) * z[0];
                
                for( Int i = 1; i < AmbDim; ++i )
                {
                    yz2 += y(k,i) * z[i];
                }
                
                value += omega[k] * std::log(std::abs( (a - two * yz2) * b ) );
            }
            
            return value * total_omega_inv;
        }
        
        
        void LineSearch_Hyperbolic_Residual()
        {
            // 2 F(0)^T.DF(0).u is the derivative of w\mapsto F(w)^T.F(w) at w = 0.
            const Real slope = two * DF.InnerProduct(F,u);
                        
            Real tau = one;
            
            const Real u_norm = u.Norm();

            // exponential map shooting from 0 to tau * u.
            Times( tau * tanhc(tau * u_norm), u, z );
            
            // Shift the point z along -w to get new updated point w .
            InverseShift(z);
            
            // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
            Shift();
            
            const Real squared_residual_at_0 = squared_residual;
            
            DifferentialAndHessian_Hyperbolic();
            
            if( linesearchQ )
            {
                
                Int backtrackings = 0;
                
                // Armijo condition for squared residual.

                ArmijoQ = squared_residual - squared_residual_at_0 - settings.Armijo_slope_factor * tau * slope < static_cast<Real>(0);
    
                while( !ArmijoQ && (backtrackings < settings.max_backtrackings) )
                {
                    ++backtrackings;
                    
                    // Estimate step size from quadratic fit if applicable.
                    
                    const Real tau_1 = settings.Armijo_shrink_factor * tau;
                    const Real tau_2 = - half * settings.Armijo_slope_factor * tau * tau * slope / ( squared_residual  - squared_residual_at_0 - tau * slope );
                    
                    tau = std::max( tau_1, tau_2 );
                                        
                    Times( tau * tanhc(tau * u_norm), u, z );
                    
                    // Shift the point z along -w to get new updated point w .
                    InverseShift(z);
                    
                    // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
                    Shift();
                    
                    DifferentialAndHessian_Hyperbolic();

                    
                    ArmijoQ = squared_residual - squared_residual_at_0 - settings.Armijo_slope_factor * tau * slope < 0;
                    
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

                const Real gamma = settings.Armijo_shrink_factor;
                
                const Real sigma = settings.Armijo_slope_factor;
                
                const Real Dphi_0 = g_factor * Dot(F,u);
                
                Int backtrackings = 0;

                // Compute potential and check Armijo condition.

//                const Real phi_0 = Potential(o);
                
                Real phi_tau = Potential(z);

                ArmijoQ = phi_tau /*- phi_0*/ - sigma * tau * Dphi_0 < 0;


                while( !ArmijoQ && (backtrackings < settings.max_backtrackings) )
                {
                    ++backtrackings;

                    const Real tau_1 = gamma * tau;
                    
                    // Estimate step size from quadratic fit if applicable.
                    const Real tau_2 = - half * sigma * tau * tau * Dphi_0 / ( phi_tau /*- phi_0*/ - tau * Dphi_0 );

                    tau = std::max( tau_1, tau_2 );
                    
                    Times( tau * tanhc(tau * u_norm), u, z );
                    
                    phi_tau = Potential(z);
                    
                    ArmijoQ = phi_tau  /*- phi_0*/ - sigma * tau * Dphi_0 < 0;
                }
            }

            // Shift the point z along -w to get new updated point w .
            InverseShift(z);

            // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
            Shift();
        }
        
        void DifferentialAndHessian_Hyperbolic()
        {
            // CAUTION: We use a different sign convention as in the paper!
            // Assemble  F = -1/2 y * omega.
            // Assemble DF = nabla F + regulatization:
            // DF_{ij} = \delta_{ij} - \sum_k x_{k,i} x_{k,j} \omega_k.

            {
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = omega[0] * y(0,i);

                    F(i) = - factor;

                    for( Int j = i; j < AmbDim; ++j )
                    {
                        DF(i,j) = - factor * y(0,j);
                    }
                }
            }
            
            for( Int k = 1; k < edge_count; ++k )
            {
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = omega[k] * y(k,i);

                    F(i) -= factor;

                    for( Int j = i; j < AmbDim; ++j )
                    {
                        DF(i,j) -= factor * y(k,j);
                    }
                }
            }
            
            // Normalize for case that the weights in omega do not sum to 1.
            for( Int i = 0; i < AmbDim; ++i )
            {
                F(i) *= total_omega_inv;

                for( Int j = i; j < AmbDim; ++j )
                {
                    DF(i,j) *= total_omega_inv;
                }
            }
            
            squared_residual = Dot(F,F);

            residual = std::sqrt( squared_residual );

            F *= half;

            // Better add the identity afterwards for precision reasons.
            for( Int i = 0; i < AmbDim; ++i )
            {
                DF(i,i) += one;
            }
        }
            
        void SearchDirection_Hyperbolic()
        {
            // Make decisions whether to continue.
            if( residual < static_cast<Real>(100.) * settings.tolerance )
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
                    continueQ = (errorestimator >settings.tolerance);
                    succeededQ = !continueQ;
                }
                else
                {
                    errorestimator = infty;
                    linesearchQ = settings.Armijo_slope_factor > zero;
                    //There is no way to reduce the residual below machine epsilon. If the algorithm reaches here, the problem is probably too ill-conditioned to be solved in machine precision.
                    continueQ = residual > settings.give_up_tolerance;
                }
            }
            else
            {
                q = big_one;
                lambda_min = eps;
                errorestimator = infty;
                linesearchQ = settings.Armijo_slope_factor > zero;
                continueQ = residual>std::max( settings.give_up_tolerance, settings.tolerance );
            }
            
            const Real c = settings.regularization * squared_residual;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                for( Int j = i; j < AmbDim; ++j )
                {
                    L(i,j) = DF(i,j) + static_cast<Real>(i==j) * c;
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
        
    public:
        
        void InverseShift( const Vector_T & s )
        {
            // Shift point w along -s.

            const Real ws2 = two * Dot(w,s);
            const Real ww  = Dot(w,w);
            const Real ss  = Dot(s,s);

            const Real a = one - ww;
            const Real b = one + ss + ws2;
            const Real c = big_one + ws2 + ww * ss;
            const Real d = one / c;

            for( Int i = 0; i < AmbDim; ++i )
            {
                w[i] = ( a * s[i] + b * w[i] ) * d;
            }
        }
        
        void Shift()
        {
            S.Shift( x.data(), w.data(), y.data(), edge_count, one );
        }
        
    public:

        
        virtual const SpherePoints_T & InitialEdgeCoordinates() const override
        {
            return x;
        }
        
        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, bool normalize = true ) override
        {
            x.Read( x_in );
            
            if( normalize )
            {
                for( Int k = 0; k < edge_count; ++k )
                {
                    Real r2 = 0;
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        r2 += x(k,i) * x(k,i);
                    }
                    
                    const Real scale = one/std::sqrt(r2);
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        x(k,i) *= scale;
                    }
                }
            }
        }
        
        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, const Int k, bool normalize = true ) override
        {
            ReadInitialEdgeCoordinates( &x_in[ AmbDim * edge_count * k], normalize);
        }
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out ) const override
        {
            x.Write(x_out);
        }
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out, const Int k ) const override
        {
            WriteInitialEdgeCoordinates( &x_out[ AmbDim * edge_count * k ]);
        }
        
        
//        virtual SpherePoints_T & EdgeCoordinates() override
//        {
//            return y;
//        }
        
        virtual const SpherePoints_T & EdgeCoordinates() const override
        {
            return y;
        }
        
        virtual void ReadEdgeCoordinates( const Real * const y_in ) override
        {
            y.Read(y_in);
        }
        
        virtual void ReadEdgeCoordinates( const Real * const y_in, const Int k ) override
        {
            ReadEdgeCoordinates( &y_in[ AmbDim * edge_count * k ]);
        }
        
        virtual void WriteEdgeCoordinates( Real * y_out ) const override
        {
            y.Write(y_out);
        }
        
        virtual void WriteEdgeCoordinates( Real * y_out, const Int k ) const override
        {
            WriteEdgeCoordinates( &y_out[ AmbDim * edge_count * k ]);
        }
        
        
        
        virtual const SpacePoints_T & SpaceCoordinates() const override
        {
            return p;
        }
                
        virtual void WriteSpaceCoordinates( Real * p_out ) const  override
        {
            p.Write(p_out);
        }
        
        virtual void WriteSpaceCoordinates( Real * p_out, const Int k ) const override
        {
            WriteSpaceCoordinates( &p_out[ (edge_count+1) * AmbDim * k ]);
        }
        
        void ComputeSpaceCoordinates() const override
        {
            //Caution: This gives only have the weight to the end vertices of the chain.
            //Thus this is only really the barycenter, if the chain is closed!
            
            Real barycenter        [AmbDim] = {};
            Real point_accumulator [AmbDim] = {};
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Real omega_k = omega[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real offset = omega_k * y(k,i);
                    
                    barycenter[i] += (point_accumulator[i] + half * offset);
                    
                    point_accumulator[i] += offset;
                }
            }

            for( Int i = 0; i < AmbDim; ++i )
            {
                p(0,i) = -barycenter[i]/edge_count;
            }

            for( Int k = 0; k < edge_count; ++k )
            {
                const Real omega_k = omega[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    p(k+1,i) = p(k,i) + omega_k * y(k,i);
                }
            }
        }
        
        
        virtual const Weights_T & Omega() const override
        {
            return omega;
        }
        
        virtual void ReadOmega( const Real * const omega_in ) override
        {
            omega.Read(omega_in);
            
            total_omega_inv = one / omega.Total();
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
            Real w_ [AmbDim] = {};
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Real weight_k = omega[k];
                for( Int i = 0; i < AmbDim; ++i )
                {
                    w_[i] += x(k,i) * weight_k;
                }
            }
            
            // Normalize in that case that omega does not sum up to 1.
            for( Int i = 0; i < AmbDim; ++i )
            {
                w_[i] *= total_omega_inv;
            }
            
            w.Read(&w_[0]);
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
            w.Write(w_out);
        }
        
        virtual void WriteShiftVector( Real * w_out, const Int k ) const override
        {
            w.Write(&w_out[ AmbDim * k]);
        }
        
        const Vector_T & ShiftVector() const
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
            return settings.max_iter;
        }
        
        
    public:
        
        virtual void OptimizeBatch(
            const Real * const x_in,
                  Real *       w_out,
                  Real *       y_out,
            const Int smaple_count,
            const Int thread_count_ = 0,
            bool normalize = true
        ) override
        {
            ptic(ClassName()+"OptimizeBatch");
            
            const Int thread_count = (thread_count_==0) ? omp_get_num_threads() : thread_count_;
            
            Tensor1<Int,Int> job_ptr = BalanceWorkLoad( smaple_count, thread_count );
            
            #pragma omp parallel num_threads( thread_count )
            {
                const Int thread = omp_get_thread_num();
                
                const Int k_begin = job_ptr[thread];
                const Int k_end   = job_ptr[thread+1];

                CLASS W( edge_count, settings );
                
                W.ReadOmega( Omega().data() );
                
                for( Int k = k_begin; k < k_end; ++k )
                {
                    W.ReadInitialEdgeCoordinates( x_in, k, normalize );

                    W.ComputeShiftVector();
                    
                    W.Optimize();
                    
                    W.WriteShiftVector( w_out, k );

                    W.WriteEdgeCoordinates( y_out, k );
                }
            }
            
            ptoc(ClassName()+"OptimizeBatch");
        }
        
        virtual void RandomClosedPolygons(
                  Real * const restrict x_out,
                  Real * const restrict w_out,
                  Real * const restrict y_out,
                  Real * const restrict K_edge_space,
                  Real * const restrict K_edge_quotient_space,
            const Int sample_count,
            const Int thread_count_ = 0
        ) const override
        {
            ptic(ClassName()+"RandomClosedPolygons");

            const Int thread_count = (thread_count_<=0) ? omp_get_num_threads() : thread_count_;

            Tensor1<Int,Int> job_ptr = BalanceWorkLoad(sample_count, thread_count);

//            valprint( "dimension   ", AmbDim       );
//            valprint( "edge_count  ", edge_count   );
//            valprint( "sample_count", sample_count );
//            valprint( "thread_count", thread_count );

            #pragma omp parallel num_threads( thread_count )
            {
                const Int thread = omp_get_thread_num();
                
                const Int k_begin = job_ptr[thread  ];
                const Int k_end   = job_ptr[thread+1];

                CLASS W( edge_count, settings );

                W.ReadOmega( Omega().data() );
                W.ReadRho( Rho().data() );

                for( Int k = k_begin; k < k_end; ++k )
                {
                    W.RandomizeInitialEdgeCoordinates();

                    W.WriteInitialEdgeCoordinates(x_out, k);
                    
                    W.ComputeShiftVector();

                    W.Optimize();
                    
                    W.WriteShiftVector(w_out, k);
                    
                    W.WriteEdgeCoordinates(y_out, k);
                    
                    W.ComputeEdgeSpaceSamplingWeight();

                    W.ComputeEdgeQuotientSpaceSamplingCorrection();

                    K_edge_space[k] = W.EdgeSpaceSamplingWeight();

                    K_edge_quotient_space[k] = W.EdgeQuotientSpaceSamplingWeight();
                }
            }
            
            ptoc(ClassName()+"::RandomClosedPolygons");
        }
        
        
        void Sample(
            Real * restrict bins_out,
            const Int bin_count_,
            Real * restrict moments_out,
            const Int moment_count_,
            Real * restrict ranges_out,
            const std::vector< std::unique_ptr<RandomVariable<AmbDim,Real,Int>> > & F_list_,
            const Int sample_count,
            const Int thread_count_ = 0
        ) const
        {
            ptic(ClassName()+"Sample");

            const Int fun_count = static_cast<Int>(F_list_.size());
            
            const Int moment_count = std::max( static_cast<Int>(3), moment_count_ );

            const Int bin_count = std::max( bin_count_, static_cast<Int>(1) );

            const Int thread_count = (thread_count_<=0) ? omp_get_num_threads() : thread_count_;

            Tensor1<Int,Int> job_ptr = BalanceWorkLoad(sample_count, thread_count);

            valprint( "dimension   ", AmbDim       );
            valprint( "edge_count  ", edge_count );
            valprint( "sample_count", sample_count );
            valprint( "fun_count   ", fun_count );
            valprint( "bin_count   ", bin_count );
            valprint( "moment_count", moment_count );
            valprint( "thread_count", thread_count );
            

            Tensor3<Real,Int> bins_global   ( 3, fun_count, bin_count,    zero );
            Tensor2<Real,Int> ranges        (    fun_count, 2                  );
            Tensor3<Real,Int> moments_global( 3, fun_count, moment_count, zero );
            Tensor1<Real,Int> factor        (    fun_count                     );

            print("Sampling the following random variables:");
            for( Int i = 0; i < fun_count; ++ i )
            {
                ranges(i,0) = F_list_[i]->MinValue(*this);
                ranges(i,1) = F_list_[i]->MaxValue(*this);

                factor(i) = static_cast<Real>(bin_count) / ( ranges(i,1) - ranges(i,0) );
                
                print("    " + F_list_[i]->Tag());
            }

            const Int lower = static_cast<Int>(0);
            const Int upper = static_cast<Int>(bin_count-1);

            #pragma omp parallel num_threads( thread_count )
            {
                const Int thread = omp_get_thread_num();

                const Int repetitions = (job_ptr[thread+1] - job_ptr[thread]);

                CLASS W( edge_count, settings );

                W.ReadOmega( Omega().data() );
                W.ReadRho( Rho().data() );

                std::vector< std::unique_ptr<RandomVariable<AmbDim,Real,Int>> > F_list;

                for( Int i = 0; i < fun_count; ++ i )
                {
                    F_list.push_back( F_list_[i]->Clone() );
                }

                Tensor3<Real,Int> bins_local   ( 3, fun_count, bin_count,    zero );
                Tensor3<Real,Int> moments_local( 3, fun_count, moment_count, zero );

                for( Int k = 0; k < repetitions; ++k )
                {
                    W.RandomizeInitialEdgeCoordinates();

                    W.ComputeShiftVector();

                    W.Optimize();

                    W.ComputeSpaceCoordinates();
                    
                    W.ComputeEdgeSpaceSamplingWeight();
                    
                    W.ComputeEdgeQuotientSpaceSamplingCorrection();
                    
                    const Real K = W.EdgeSpaceSamplingWeight();

                    const Real K_quot = W.EdgeQuotientSpaceSamplingWeight();

                    for( Int i = 0; i < fun_count; ++i )
                    {
                        auto & F = *F_list[i];

                        const Real val = F(W);
                        
                        Real values [3] = {static_cast<Real>(1),K,K_quot};
                        
//                        const Int bin_idx = std::clamp(
//                           static_cast<Int>(std::floor( factor[i] * (val - ranges(i,0)) )),
//                           lower,
//                           upper
//                        );
                        
                        const Int bin_idx = static_cast<Int>(std::floor( factor[i] * (val - ranges(i,0)) ));
                        
                        if( (bin_idx <= upper) && (bin_idx >= lower) )
                        {
                            bins_local(0,i,bin_idx) += static_cast<Real>(1);
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

                #pragma omp critical
                {
//                    valprint("thread",thread);
//                    print( "data_local = " + data_local.ToString() );
                    for( Int l = 0; l < 3; ++l )
                    {
                        for( Int i = 0; i < fun_count; ++i )
                        {
                            for( Int j = 0; j < bin_count; ++j )
                            {
                                bins_global(l,i,j) += bins_local(l,i,j);
                            }
                            
                            for( Int j = 0; j < moment_count; ++j )
                            {
                                moments_global(l,i,j) += moments_local(l,i,j);
                            }
                        }
                    }
                }
            }

            bins_global.Write( bins_out );
            moments_global.Write( moments_out );
            ranges.Write( ranges_out );
            
            ptoc(ClassName()+"::Sample");
        }
        
        virtual void Sample(
            Real * restrict bins_out,
            const Int bin_count_,
            Real * restrict moments_out,
            const Int moment_count_,
            Real * restrict ranges_out,
            const std::vector< std::unique_ptr<RandomVariableBase<Real,Int>> > & F_list_,
            const Int sample_count,
            const Int thread_count_ = 0
        ) const override
        {
            ptic(ClassName()+"Sample (polymorphic)");

            const Int fun_count = static_cast<Int>(F_list_.size());
            
            const Int moment_count = std::max( static_cast<Int>(3), moment_count_ );

            const Int bin_count = std::max( bin_count_, static_cast<Int>(1) );

            const Int thread_count = (thread_count_<=0) ? omp_get_num_threads() : thread_count_;

            Tensor1<Int,Int> job_ptr = BalanceWorkLoad(sample_count, thread_count);

            valprint( "dimension   ", AmbDim       );
            valprint( "edge_count  ", edge_count );
            valprint( "sample_count", sample_count );
            valprint( "fun_count   ", fun_count );
            valprint( "bin_count   ", bin_count );
            valprint( "moment_count", moment_count );
            valprint( "thread_count", thread_count );
            

            Tensor3<Real,Int> bins_global   ( 3, fun_count, bin_count,    zero );
            Tensor2<Real,Int> ranges        (    fun_count, 2                  );
            Tensor3<Real,Int> moments_global( 3, fun_count, moment_count, zero );
            Tensor1<Real,Int> factor        (    fun_count                     );

            print("Sampling the following random variables:");
            for( Int i = 0; i < fun_count; ++ i )
            {
                const size_t i_ = static_cast<size_t>(i);
                ranges(i,0) = F_list_[i_]->MinValue(*this);
                ranges(i,1) = F_list_[i_]->MaxValue(*this);

                factor(i) = static_cast<Real>(bin_count) / ( ranges(i,1) - ranges(i,0) );
                
                print("    " + F_list_[i_]->Tag());
            }

            const Int lower = static_cast<Int>(0);
            const Int upper = static_cast<Int>(bin_count-1);

            #pragma omp parallel num_threads( thread_count )
            {
                const Int thread = omp_get_thread_num();

                const Int repetitions = (job_ptr[thread+1] - job_ptr[thread]);

                CLASS W( edge_count, settings );

                W.ReadOmega( Omega().data() );
                W.ReadRho( Rho().data() );

                std::vector< std::unique_ptr<RandomVariableBase<Real,Int>> > F_list;

                for( Int i = 0; i < fun_count; ++ i )
                {
                    F_list.push_back( F_list_[static_cast<size_t>(i)]->Clone() );
                }

                Tensor3<Real,Int> bins_local   ( 3, fun_count, bin_count,    zero );
                Tensor3<Real,Int> moments_local( 3, fun_count, moment_count, zero );

                for( Int k = 0; k < repetitions; ++k )
                {
                    W.RandomizeInitialEdgeCoordinates();

                    W.ComputeShiftVector();

                    W.Optimize();

                    W.ComputeSpaceCoordinates();
                    
                    W.ComputeEdgeSpaceSamplingWeight();
                    
                    W.ComputeEdgeQuotientSpaceSamplingCorrection();
                    
                    const Real K = W.EdgeSpaceSamplingWeight();

                    const Real K_quot = W.EdgeQuotientSpaceSamplingWeight();

                    for( Int i = 0; i < fun_count; ++i )
                    {
                        auto & F = *F_list[static_cast<size_t>(i)];

                        const Real val = F(W);
                        
                        Real values [3] = {static_cast<Real>(1),K,K_quot};
                        
//                        const Int bin_idx = std::clamp(
//                           static_cast<Int>(std::floor( factor[i] * (val - ranges(i,0)) )),
//                           lower,
//                           upper
//                        );
                        
                        const Int bin_idx = static_cast<Int>(std::floor( factor[i] * (val - ranges(i,0)) ));
                        
                        if( (bin_idx <= upper) && (bin_idx >= lower) )
                        {
                            bins_local(0,i,bin_idx) += static_cast<Real>(1);
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

                #pragma omp critical
                {
//                    valprint("thread",thread);
//                    print( "data_local = " + data_local.ToString() );
                    for( Int l = 0; l < 3; ++l )
                    {
                        for( Int i = 0; i < fun_count; ++i )
                        {
                            for( Int j = 0; j < bin_count; ++j )
                            {
                                bins_global(l,i,j) += bins_local(l,i,j);
                            }
                            
                            for( Int j = 0; j < moment_count; ++j )
                            {
                                moments_global(l,i,j) += moments_local(l,i,j);
                            }
                        }
                    }
                }
            }

            bins_global.Write( bins_out );
            moments_global.Write( moments_out );
            ranges.Write( ranges_out );
            
            ptoc(ClassName()+"::Sample (polymorphic)");
        }
        
#if defined(PLCTOPOLOGY_H)
        
        std::map<std::string, std::tuple<Real,Real,Real>> SampleHOMFLY(
            const Int sample_count,
            const Int thread_count_ = 0
        ) const
        {
            ptic(ClassName()+"SampleHOMFLY");
            
            std::map<std::string, std::tuple<Real,Real,Real>> map_global;
            
            if( AmbDim != 3 )
            {
                eprint("SampleHOMFLY is only available in 3D.");
                return map_global;
            }

            const Int thread_count = (thread_count_<=0) ? omp_get_num_threads() : thread_count_;

            Tensor1<Int,Int> job_ptr = BalanceWorkLoad(sample_count, thread_count);
            
            valprint( "dimension   ", AmbDim       );
            valprint( "edge_count  ", edge_count   );
            valprint( "sample_count", sample_count );
            valprint( "thread_count", thread_count );
            
            gsl_rng_env_setup();
            
            Tensor1<unsigned long,Int> seeds ( thread_count );
            
            std::random_device r;
            
            for( Int thread = 0 ; thread < thread_count; ++ thread )
            {
                seeds[thread] = (std::numeric_limits<unsigned int>::max()+1) * r() + r();
            }
            
            #pragma omp parallel num_threads( thread_count )
            {
                const Int thread = omp_get_thread_num();

                const Int repetitions = (job_ptr[thread+1] - job_ptr[thread]);

                CLASS W( edge_count, settings );

                W.ReadOmega( Omega().data() );
                W.ReadRho( Rho().data() );

                std::map<std::string, std::tuple<Real,Real,Real>> map_loc;

                bool openQ = false;

                int color_count = 0;

                int n = static_cast<int>(edge_count);

                gsl_rng *rng = gsl_rng_alloc( gsl_rng_mt19937 );

                gsl_rng_set( rng, seeds[thread] );
                
                plCurve * Gamma = plc_new( 1, &n, &openQ, &color_count );
//
//                auto * p = &Gamma->cp[0].vt[0];

                auto & p = W.SpaceCoordinates();
                
                for( Int l = 0; l < repetitions; ++l )
                {
//                    valprint("l",l);
                    
                    W.RandomizeInitialEdgeCoordinates();

                    W.ComputeShiftVector();

                    W.Optimize();

                    W.ComputeSpaceCoordinates();

                    W.ComputeEdgeSpaceSamplingWeight();

                    W.ComputeEdgeQuotientSpaceSamplingCorrection();
                    
                    const Real K = W.EdgeSpaceSamplingWeight();

                    const Real K_quot = W.EdgeQuotientSpaceSamplingWeight();
                    
                    auto * comp = &Gamma->cp[0];
                    
                    for( int k = 0; k < edge_count; ++k )
                    {
                        auto * v = &comp->vt[k].c[0];

                        for( int i = 0; i < 3; ++i )
                        {
                            v[i] = p(k,i);
                        }
                    }
                    
                    plc_fix_wrap( Gamma );
                    
                    
//                    char * polynomial = (char *) malloc( sizeof(char) * 3);
//                    
//                    polynomial[0] = 'a';
//                    polynomial[1] = 'b';
//                    polynomial[2] = 'c';
                    
                    char * polynomial = plc_homfly( rng, Gamma );

                    std::string s("");
                    
                    if( polynomial != nullptr )
                    {
                        s.append(polynomial);
                        free( polynomial );
                    }
                    else
                    {
                        s.append("FAILED");
                    }
                    
//                    s << polynomial;
//
//                    std::string str = s.str();
//                    print(s);
                    
                    map_loc.try_emplace(s, std::tie(zero,zero,zero) );
                    
                    auto & tuple = map_loc[s];
                    
                    std::get<0>(tuple) += one;
                    std::get<1>(tuple) += K;
                    std::get<2>(tuple) += K_quot;
                    

                }
                
                #pragma omp critical
                {
                    for ( auto const & [key, val] : map_loc )
                    {
                        map_global.try_emplace( key, std::tie(zero,zero,zero) );
                        
                        auto & from = val;
                        auto & to   = map_global[key];
                        
                        std::get<0>(to) += std::get<0>(from);
                        std::get<1>(to) += std::get<1>(from);
                        std::get<2>(to) += std::get<2>(from);
                    }
                }
                
                gsl_rng_free( rng );
                
                plc_free( Gamma );
            }


            ptoc(ClassName()+"::SampleHOMFLY");
            
            return map_global;
        }
        
#endif
        
    public:
        
        void RandomizeInitialEdgeCoordinates() override
        {
            for( Int k = 0; k < edge_count; ++k )
            {
                Real r2 = static_cast<Real>(0);

                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real z = normal_dist( random_engine );
                    
                    x(k,i) = z;
                    
                    r2 += z * z;
                }

                Real r_inv = one/std::sqrt(r2);

                for( Int i = 0; i < AmbDim; ++i )
                {
                    x(k,i) *= r_inv;
                }
            }
        }
        
        
        void RandomSphericalPoints(
            Real * restrict x_out,
            const Int sample_count,
            const Int thread_count_ = 1
        ) const override
        {
            const Int thread_count = (thread_count_<=0) ? omp_get_num_threads() : thread_count_;

            if( thread_count == 1 )
            {
                for( Int k = 0; k < edge_count; ++k )
                {
                    Real r2 = static_cast<Real>(0);

                    Real v [AmbDim];
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        v[i] = normal_dist( random_engine );

                        r2 += v[i] * v[i];
                    }

                    Real r_inv = one/std::sqrt(r2);

                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        x_out[ AmbDim * k + i ] = v[i] * r_inv;
                    }
                }
            }
            else
            {
                Tensor1<Int,Int> job_ptr = BalanceWorkLoad(sample_count, thread_count);

                #pragma omp parallel num_threads( thread_count )
                {
                    Real v [AmbDim];

                    const Int thread = omp_get_thread_num();

                    std::random_device r;
                    
                    std::seed_seq seed { r(), r(), r(), r() };
                    
                    std::mt19937_64 random_engine_loc { seed };

                    std::normal_distribution<Real> dist { static_cast<Real>(0),static_cast<Real>(1) };

                    const Int l_begin = job_ptr[thread];
                    const Int l_end   = job_ptr[thread+1];

                    for( Int l = l_begin; l < l_end; ++l )
                    {
                        Real * restrict x_ = &x_out[AmbDim * edge_count * l];

                        for( Int k = 0; k < edge_count; ++k )
                        {
                            Real r2 = static_cast<Real>(0);

                            for( Int i = 0; i < AmbDim; ++i )
                            {
                                v[i] = dist( random_engine_loc );

                                r2 += v[i] * v[i];
                            }

                            Real r_inv = one/std::sqrt(r2);

                            for( Int i = 0; i < AmbDim; ++i )
                            {
                                x_[ AmbDim * k + i ] = v[i] * r_inv;
                            }
                        }
                    }
                }
            }
        }
        
    protected:
        
        ShiftMap<AmbDim,Real,Int> & Shifter()
        {
            return S;
        }
        
    protected:
        
        Real tanhc( const Real t ) const
        {
            // Computes tanh(t)/t in a stable way by using a Padé approximation around t = 0.
            constexpr Real a0 = static_cast<Real>(1);
            constexpr Real a1 = static_cast<Real>(7)/static_cast<Real>(51);
            constexpr Real a2 = static_cast<Real>(1)/static_cast<Real>(255);
            constexpr Real a3 = static_cast<Real>(2)/static_cast<Real>(69615);
            constexpr Real a4 = static_cast<Real>(1)/static_cast<Real>(34459425);
            
            constexpr Real b0 = static_cast<Real>(1);
            constexpr Real b1 = static_cast<Real>(8)/static_cast<Real>(17);
            constexpr Real b2 = static_cast<Real>(7)/static_cast<Real>(255);
            constexpr Real b3 = static_cast<Real>(4)/static_cast<Real>(9945);
            constexpr Real b4 = static_cast<Real>(1)/static_cast<Real>(765765);
            
            const Real t2 = t * t;
            
            const Real result = ( t2 <= one ) ? (
                a0 + t2 * (a1 + t2 * (a2 + t2 * (a3 + t2 * a4)))
            )/(
                b0 + t2 * (b1 + t2 * (b2 + t2 * (b3 + t2 * b4)))
            )
            : ( t2 <= static_cast<Real>(7) ) ? std::tanh(t)/t : one/std::abs(t);
            
            return result;
        }
        
    public:
        
        Int AmbientDimension() const override
        {
            return AmbDim;
        }
        
        virtual std::string ClassName() const override
        {
            return TO_STD_STRING(CLASS)+"<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+">";
        }
    };
    
    
    template<typename Real = double, typename Int = long long>
    std::unique_ptr<BASE> MakeCyclicSampler(
            const int amb_dim,
            const Int edge_count_,
            const CyclicSamplerSettings<Real,Int> settings_ = CyclicSamplerSettings<Real,Int>()
    )
    {
        CyclicSamplerSettings<Real,Int> settings (settings_);
        switch( amb_dim )
        {
            case 2:
            {
                return std::make_unique<CLASS<2,Real,Int>>(edge_count_,settings);
            }
            case 3:
            {
                return std::make_unique<CLASS<3,Real,Int>>(edge_count_,settings);
            }
            case 4:
            {
                return std::make_unique<CLASS<4,Real,Int>>(edge_count_,settings);
            }
                
            default:
            {
                eprint("Make"+TO_STD_STRING(CLASS)+": ambient dimension "+ToString(amb_dim)+" not supported. Using default dimension 3.");
                return std::make_unique<CLASS<3,Real,Int>>(edge_count_,settings);
            }
        }
    }
        
#undef BASE
#undef CLASS
    
} // namespace CyclicSampler
