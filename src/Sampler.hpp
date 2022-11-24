#pragma once

#define x(k,i) x_buffer[AmbDim*k+i]
#define y(k,i) y_buffer[AmbDim*k+i]
#define p(k,i) p_buffer[AmbDim*k+i]

namespace CycleSampler
{
    
    template<int AmbDim, typename Real = double, typename Int = long long, bool copy = true >
    class Sampler : public SamplerBase<Real,Int>
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using Base_T            = SamplerBase<Real,Int>;
        using Vector_T          = SmallVector<AmbDim,Real,Int>;
        using SymmetricMatrix_T = SmallSymmetricMatrix<AmbDim,Real,Int>;
        
        using HyperbolicPoint_T = Vector_T;
        
        using SpherePoints_T = typename Base_T::SpherePoints_T;
        using SpacePoints_T  = typename Base_T::SpacePoints_T;
        using Weights_T      = typename Base_T::Weights_T;
        using Setting_T      = typename Base_T::Setting_T;
        
        
        using Base_T::settings;
        using Base_T::edge_count;
        using Base_T::Settings;
        using Base_T::random_engine;
        using Base_T::normal_dist;
        
        Sampler() : Base_T() {}
        
        virtual ~Sampler()
        {
            if constexpr ( copy )
            {
                safe_free( x_buffer );
                safe_free( y_buffer );
                safe_free( p_buffer );
            }
            
            safe_free( r );
            safe_free( rho );
        }
        
        explicit Sampler(
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   Base_T( edge_count_, settings_ )
        {
            if constexpr ( copy )
            {
                safe_alloc( x_buffer,  edge_count    * AmbDim );
                safe_alloc( y_buffer,  edge_count    * AmbDim );
                safe_alloc( p_buffer, (edge_count+1) * AmbDim );
            }
            
            safe_alloc ( r,   edge_count );
            fill_buffer( r,   edge_count, one / edge_count );
            safe_alloc ( rho, edge_count );
            fill_buffer( rho, edge_count, one );
            
            total_r_inv = one;
        }
        
        explicit Sampler(
            const Real * restrict const r_in,
            const Real * restrict const rho_in,
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   Base_T( edge_count_, settings_ )
        {
            if constexpr ( copy )
            {
                safe_alloc( x_buffer,  edge_count    * AmbDim );
                safe_alloc( y_buffer,  edge_count    * AmbDim );
                safe_alloc( p_buffer, (edge_count+1) * AmbDim );
            }
            
            safe_alloc ( r, edge_count );
            copy_buffer( r_in, r, edge_count );
            
            safe_alloc ( rho, edge_count );
            copy_buffer( rho_in, rho, edge_count );
        }


    protected:
        
        Real * restrict x_buffer = nullptr;
        Real * restrict y_buffer = nullptr;
        Real * restrict p_buffer = nullptr;
        
        Real * restrict r   = nullptr;
        Real * restrict rho = nullptr;
        
        Real total_r_inv = one;
        
        Vector_T F;           // right hand side of Newton iteration.
        SymmetricMatrix_T DF; // nabla F(0) with respect to measure ys

        Real w [AmbDim];           // current point in hyperbolic space.
        Real u [AmbDim];           // update direction
        Real z [AmbDim];           // Multiple purpose buffer.
        Real A [AmbDim][AmbDim];   // For Cholesky factorization.
        
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
                x_buffer, &w[0], y_buffer, r, rho, edge_count
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
        
        Real Potential()
        {
            Real value = 0;
            
            Real zz = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                zz += z[i] * z[i];
            }
            
            
            const Real a = big_one + zz;
            const Real c = (big_one-zz);
            
            const Real b = one/c;
            
            for( Int k = 0; k < edge_count; ++k )
            {
                Real yz = y(k,0) * z[0];
                
                for( Int i = 1; i < AmbDim; ++i )
                {
                    yz += y(k,i) * z[i];
                }
                
                value += r[k] * std::log(std::abs( (a - two * yz) * b ) );
            }
            
            return value * total_r_inv;
        }
        
        
        void LineSearch_Hyperbolic_Residual()
        {
            // 2 F(0)^T.DF(0).u is the derivative of w\mapsto F(w)^T.F(w) at w = 0.
            const Real slope = two * DF.InnerProduct(F,u);
                        
            Real tau = one;
            
            Real uu = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                uu += u[i] * u[i];
            }
            
            const Real u_norm = std::sqrt(uu);
            

            // exponential map shooting from 0 to tau * u.
            {
                const Real scale = tau * tanhc(tau * u_norm);
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    z[i] = scale * u[i];
                }
            }
            
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

            Real uu = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                uu += u[i] * u[i];
            }
            
            const Real u_norm = std::sqrt(uu);

            // exponential map shooting from 0 to tau * u.
            {
                const Real scale = tau * tanhc(tau * u_norm);
                for( Int i = 0; i < AmbDim; ++i )
                {
                    z[i] = scale * u[i];
                }
            }
            
            if( linesearchQ )
            {
                //Linesearch with potential as merit function.

                const Real gamma = settings.Armijo_shrink_factor;
                
                const Real sigma = settings.Armijo_slope_factor;
                
                Real Dphi_0 = 0;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    Dphi_0 += F[i] * u[i];
                }
                
                Dphi_0 *= g_factor;
                
                Int backtrackings = 0;

                // Compute potential and check Armijo condition.

//                const Real phi_0 = Potential(o);
                
                Real phi_tau = Potential();

                ArmijoQ = phi_tau /*- phi_0*/ - sigma * tau * Dphi_0 < 0;


                while( !ArmijoQ && (backtrackings < settings.max_backtrackings) )
                {
                    ++backtrackings;

                    const Real tau_1 = gamma * tau;
                    
                    // Estimate step size from quadratic fit if applicable.
                    const Real tau_2 = - half * sigma * tau * tau * Dphi_0 / ( phi_tau /*- phi_0*/ - tau * Dphi_0 );

                    tau = std::max( tau_1, tau_2 );
                    
                    {
                        const Real scale = tau * tanhc(tau * u_norm);
                        for( Int i = 0; i < AmbDim; ++i )
                        {
                            z[i] = scale * u[i];
                        }
                    }
                    
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

            {
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = r[0] * y(0,i);

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
                    const Real factor = r[k] * y(k,i);

                    F(i) -= factor;

                    for( Int j = i; j < AmbDim; ++j )
                    {
                        DF(i,j) -= factor * y(k,j);
                    }
                }
            }
            
            // Normalize for case that the weights in r do not sum to 1.
            for( Int i = 0; i < AmbDim; ++i )
            {
                F(i) *= total_r_inv;

                for( Int j = i; j < AmbDim; ++j )
                {
                    DF(i,j) *= total_r_inv;
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
                    A[i][j] = DF(i,j) + static_cast<Real>(i==j) * c;
                }
            }
            
            Cholesky();
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                u[i] = -F[i];
            }
            
            CholeskySolve();
        }
        
        void Gradient_Hyperbolic()
        {
            for( Int i = 0; i < AmbDim; ++i )
            {
                u[i] = -g_factor_inv * F[i];
            }
        }
        
        void Gradient_Planar()
        {
            // The factor 2 is here to reproduce the Abikoff-Ye algorithm (in the absence of linesearch.)
            for( Int i = 0; i < AmbDim; ++i )
            {
                u[i] = - F[i];
            }
        }
        
    public:
        
        void InverseShift()
        {
            // Shift point w along -z.

            Real ww  = 0;
            Real wz2 = 0;
            Real zz  = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                ww  += w[i] * w[i];
                wz2 += w[i] * z[i];
                zz  += z[i] * z[i];
            }
            
            wz2 *= two;

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
            S.Shift( x_buffer, &w[0], y_buffer, edge_count, one );
        }
        
    public:

        void Normalize()
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
        
        virtual const Real * InitialEdgeCoordinates() const override
        {
            return x_buffer;
        }
        
        virtual void LoadInitialEdgeCoordinates( Real * const x_in, bool normalize = true ) override
        {
            if constexpr ( copy )
            {
                copy_buffer( x_in, x_buffer, edge_count * AmbDim );
            }
            else
            {
                x_buffer = x_in;
            }
            
            if( normalize )
            {
                Normalize();
            }
        }
        
        virtual void LoadInitialEdgeCoordinates( Real * const x_in, const Int k, bool normalize = true ) override
        {
            LoadInitialEdgeCoordinates( &x_in[ AmbDim * edge_count * k], normalize);
        }
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out ) const override
        {
            copy_buffer( x_buffer, x_out, edge_count * AmbDim );
        }
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out, const Int k ) const override
        {
            WriteInitialEdgeCoordinates( &x_out[ AmbDim * edge_count * k ]);
        }
        
        
//        virtual SpherePoints_T & EdgeCoordinates() override
//        {
//            return y;
//        }
        
        virtual const Real * EdgeCoordinates() const override
        {
            return y_buffer;
        }
        
        virtual void LoadEdgeCoordinates( Real * const y_in ) override
        {
            if constexpr ( copy )
            {
                copy_buffer( y_in, y_buffer, edge_count * AmbDim );
            }
            else
            {
                y_buffer = y_in;
            }
        }
        
        virtual void LoadEdgeCoordinates( Real * const y_in, const Int k ) override
        {
            LoadEdgeCoordinates( &y_in[ AmbDim * edge_count * k ]);
        }
        
        virtual void WriteEdgeCoordinates( Real * y_out ) const override
        {
            copy_buffer( y_buffer, y_out, edge_count * AmbDim );
        }
        
        virtual void WriteEdgeCoordinates( Real * y_out, const Int k ) const override
        {
            WriteEdgeCoordinates( &y_out[ AmbDim * edge_count * k ]);
        }
        
        
        
        virtual const Real * SpaceCoordinates() const override
        {
            return p_buffer;
        }
                
        virtual void WriteSpaceCoordinates( Real * p_out ) const  override
        {
            copy_buffer( p_buffer, p_out, (edge_count+1) * AmbDim );
        }
        
        virtual void WriteSpaceCoordinates( Real * p_out, const Int k ) const override
        {
            WriteSpaceCoordinates( &p_out[ (edge_count+1) * AmbDim * k ]);
        }
        
        virtual void ComputeSpaceCoordinates() const override
        {
            //Caution: This gives only have the weight to the end vertices of the chain.
            //Thus this is only really the barycenter, if the chain is closed!
            
            Real barycenter        [AmbDim] = {};
            Real point_accumulator [AmbDim] = {};
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Real r_k = r[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real offset = r_k * y(k,i);
                    
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
                const Real r_k = r[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    p(k+1,i) = p(k,i) + r_k * y(k,i);
                }
            }
        }
        

        virtual const Real * EdgeLengths() const override
        {
            return r;
        }
        
        virtual void ReadEdgeLengths( const Real * const r_in ) override
        {
            copy_buffer( r_in, r, edge_count );
            
            Real sum = 0;
            
            for( Int i = 0; i < edge_count; ++i )
            {
                sum += r[i];
            }
            
            total_r_inv = one / sum;
        }
        
        
        virtual const Real * Rho() const override
        {
            return rho;
        }
        
        virtual void ReadRho( const Real * const rho_in ) override
        {
            copy_buffer( rho_in, rho, edge_count );
        }
        
        
        virtual void ComputeShiftVector() override
        {
            Real w_ [AmbDim] = {};
            
            for( Int k = 0; k < edge_count; ++k )
            {
                const Real r_k = r[k];
                for( Int i = 0; i < AmbDim; ++i )
                {
                    w_[i] += x(k,i) * r_k;
                }
            }
            
            // Normalize in that case that r does not sum up to 1.
            for( Int i = 0; i < AmbDim; ++i )
            {
                w[i] = w_[i] * total_r_inv;
            }
        }
        
        virtual void ReadShiftVector( const Real * const w_in ) override
        {
            Real ww = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                w[i] = w_in[i];
                
                ww += w[i] * w[i];
            }
            
            // Use Euclidean barycenter as initial guess if the supplied initial guess does not make sense.
            if( ww > small_one )
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
            for( Int i = 0; i < AmbDim; ++i )
            {
                w_out[i] = w[i];
            }
        }
        
        virtual void WriteShiftVector( Real * w_out, const Int k ) const override
        {
            WriteShiftVector(&w_out[ AmbDim * k]);
        }
        
        const Real * ShiftVector() const
        {
            return &w[0];
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
                  Real * const x_in,
                  Real * const w_out,
                  Real * const y_out,
            const Int sample_count,
            const Int thread_count = 1,
            bool normalize = true
        ) override
        {
            ptic(ClassName()+"OptimizeBatch");
            
            JobPointers<Int> job_ptr ( sample_count, thread_count );
            
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                const Int k_begin = job_ptr[thread];
                const Int k_end   = job_ptr[thread+1];

                Sampler W( edge_count, settings );
                
                W.ReadEdgeLengths( EdgeLengths() );
                
                for( Int k = k_begin; k < k_end; ++k )
                {
                    W.LoadInitialEdgeCoordinates( x_in, k, normalize );
                    
                    if constexpr ( !copy )
                    {
                        W.LoadEdgeCoordinates( y_out, k );
                    }
                    
                    W.ComputeShiftVector();
                    
                    W.Optimize();
                    
                    W.WriteShiftVector( w_out, k );
                    
                    if constexpr ( copy )
                    {
                        W.WriteEdgeCoordinates( y_out, k );
                    }
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
            const Int thread_count = 1
        ) const override
        {
            ptic(ClassName()+"RandomClosedPolygons");

            JobPointers<Int> job_ptr ( sample_count, thread_count );

//            valprint( "dimension   ", AmbDim       );
//            valprint( "edge_count  ", edge_count   );
//            valprint( "sample_count", sample_count );
//            valprint( "thread_count", thread_count );

            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                const Int k_begin = job_ptr[thread  ];
                const Int k_end   = job_ptr[thread+1];

                Sampler W( edge_count, settings );

                W.ReadEdgeLengths( EdgeLengths() );
                W.ReadRho( Rho() );

                for( Int k = k_begin; k < k_end; ++k )
                {
                    if constexpr ( !copy )
                    {
                        W.LoadInitialEdgeCoordinates( x_out, k );
                        W.LoadEdgeCoordinates( y_out, k );
                    }
                    
                    W.RandomizeInitialEdgeCoordinates();

                    if constexpr ( copy )
                    {
                        W.WriteInitialEdgeCoordinates( x_out, k );
                    }
                    
                    W.ComputeShiftVector();

                    W.Optimize();
                    
                    W.WriteShiftVector( w_out, k );
                    
                    if constexpr ( copy )
                    {
                        W.WriteEdgeCoordinates(y_out, k);
                    }
                    
                    W.ComputeEdgeSpaceSamplingWeight();

                    W.ComputeEdgeQuotientSpaceSamplingCorrection();

                    K_edge_space[k] = W.EdgeSpaceSamplingWeight();

                    K_edge_quotient_space[k] = W.EdgeQuotientSpaceSamplingWeight();
                }
            }
            
            ptoc(ClassName()+"::RandomClosedPolygons");
        }
        
        
        
        
        // moments: A 3D-array of size 3 x fun_count x bin_count. Entry moments(i,j,k) will store the sampled weighted k-th moment of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i (see above).
        // ranges: Specify the range for binning: For j-th function in F_list, the range from ranges(j,0) to ranges(j,1) will be devided into bin_count bins. The user is supposed to provide meaningful ranges. Some rough guess might be obtained by calling the random variables on the prepared CyclicSampler_T C.
        virtual void Sample_Binned(
            Real * restrict bins_out,
            const Int bin_count_,
            Real * restrict moments_out,
            const Int moment_count_,
            const Real * restrict ranges,
            const std::vector< std::unique_ptr<RandomVariableBase<Real,Int>> > & F_list_,
            const Int sample_count,
            const Int thread_count = 1
        ) const override
        {
            ptic(ClassName()+"Sample_Binned (polymorphic)");

            const Int fun_count = static_cast<Int>(F_list_.size());
            
            const Int moment_count = std::max( static_cast<Int>(3), moment_count_ );

            const Int bin_count = std::max( bin_count_, static_cast<Int>(1) );

            JobPointers<Int> job_ptr ( sample_count, thread_count );

            valprint( "dimension   ", AmbDim       );
            valprint( "edge_count  ", edge_count );
            valprint( "sample_count", sample_count );
            valprint( "fun_count   ", fun_count );
            valprint( "bin_count   ", bin_count );
            valprint( "moment_count", moment_count );
            valprint( "thread_count", thread_count );
            

            Tensor3<Real,Int> bins_global   ( bins_out,    3, fun_count, bin_count    );
            Tensor3<Real,Int> moments_global( moments_out, 3, fun_count, moment_count );
            Tensor1<Real,Int> factor        (                 fun_count               );
            
            print("Sampling the following random variables:");
            for( Int i = 0; i < fun_count; ++ i )
            {
                const size_t i_ = static_cast<size_t>(i);
                factor(i) = static_cast<Real>(bin_count) / ( ranges[2*i+1] - ranges[2*i+0] );
                
                print("    " + F_list_[i_]->Tag());
            }

            const Int lower = static_cast<Int>(0);
            const Int upper = static_cast<Int>(bin_count-1);

            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                const Int repetitions = (job_ptr[thread+1] - job_ptr[thread]);

                Sampler W( edge_count, settings );

                W.ReadEdgeLengths( EdgeLengths() );
                W.ReadRho( Rho() );

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
                        
                        const Int bin_idx = static_cast<Int>(std::floor( factor[i] * (val - ranges[2*i]) ));
                        
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
            
            ptoc(ClassName()+"::Sample_Binned (polymorphic)");
        }
        
        virtual void NormalizeBinnedSamples(
            Real * restrict bins,
            const Int bin_count,
            Real * restrict moments,
            const Int moment_count,
            const Int fun_count
        ) const override
        {
            ptic(ClassName()+"::NormalizeBinnedSamples");
            for( Int i = 0; i < 3; ++i )
            {
                for( Int j = 0; j < fun_count; ++j )
                {
                    // Normalize bins and moments.
                    
                    Real * restrict const bins_i_j = &bins[ (i*fun_count+j)*bin_count ];
                    
                    Real * restrict const moments_i_j = &moments[ (i*fun_count+j)*moment_count ];
                    
                    // The field for zeroth moment is assumed to contain the total mass.
                    Real factor = Real(1)/moments_i_j[0];

                    
                    for( Int k = 0; k < bin_count; ++k )
                    {
                        bins_i_j[k] *= factor;
                    }
               
                    for( Int k = 0; k < moment_count; ++k )
                    {
                        moments_i_j[k] *= factor;
                    }
                }
            }
            ptoc(ClassName()+"::NormalizeBinnedSamples");
        }
        
#if defined(PLCTOPOLOGY_H)
        
        std::map<std::string, std::tuple<Real,Real,Real>> SampleHOMFLY(
            const Int sample_count,
            const Int thread_count = 1
        ) const
        {
            ptic(ClassName()+"SampleHOMFLY");
            
            std::map<std::string, std::tuple<Real,Real,Real>> map_global;
            
            if( AmbDim != 3 )
            {
                eprint("SampleHOMFLY is only available in 3D.");
                return map_global;
            }

            JobPointers<Int> job_ptr (sample_count, thread_count);
            
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
            
            #pragma omp parallel for num_threads( thread_count )
            for( Int thread = 0; thread < thread_count; ++thread )
            {
                const Int repetitions = (job_ptr[thread+1] - job_ptr[thread]);

                Sampler W( edge_count, settings );

                W.ReadEdgeLengths( EdgeLengths().data() );
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
                JobPointers<Int> job_ptr ( sample_count, thread_count );

                #pragma omp parallel for num_threads( thread_count )
                for( Int thread = 0; thread < thread_count; ++thread )
                {
                    Real v [AmbDim];

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
        
        void Cholesky()
        {
            for( Int k = 0; k < AmbDim; ++k )
            {
                const Real a = A[k][k] = std::sqrt(A[k][k]);
                const Real ainv = one/a;

                for( Int j = k+1; j < AmbDim; ++j )
                {
                    A[k][j] *= ainv;
                }

                for( Int i = k+1; i < AmbDim; ++i )
                {
                    for( Int j = i; j < AmbDim; ++j )
                    {
                        A[i][j] -= A[k][i] * A[k][j];
                    }
                }
            }
        }
        
        void CholeskySolve()
        {
            //In-place solve.
            
            // Lower triangular back substitution
            for( Int i = 0; i < AmbDim; ++i )
            {
                for( Int j = 0; j < i; ++j )
                {
                    u[i] -= A[j][i] * u[j];
                }
                u[i] /= A[i][i];
            }
            
            // Upper triangular back substitution
            for( Int i = AmbDim-1; i > -1; --i )
            {
                for( Int j = i+1; j < AmbDim; ++j )
                {
                    u[i] -= A[i][j] * u[j];
                }
                u[i] /= A[i][i];
            }
        }
        
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
            return "Sampler<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+">";
        }
    };
    
} // namespace CycleSampler
  

#undef x
#undef y
#undef p


//#pragma once
//
//
//#define x(k,i) x_buffer[AmbDim*k+i]
//#define y(k,i) y_buffer[AmbDim*k+i]
//#define p(k,i) p_buffer[AmbDim*k+i]
//
//namespace CycleSampler
//{
//
//    template<int AmbDim, typename Real = double, typename Int = long long>
//    class Sampler : public SamplerBase<Real,Int>
//
//    {
//        ASSERT_FLOAT(Real);
//        ASSERT_INT(Int);
//
//    public:
//
//        using Base_T = SamplerBase<Real,Int>;
//
////        using SpherePoints_T = typename Base_T::SpherePoints_T;
////        using SpacePoints_T  = typename Base_T::SpacePoints_T;
//        using Setting_T      = typename Base_T::Setting_T;
//
//        using RandomVariableBase_T = RandomVariableBase<Real,Int>;
//        using RandomVariable_T     = RandomVariable<AmbDim,Real,Int>;
//
//        using Vector_T  = SmallVector<AmbDim, Real, Int>;
//
//        using Base_T::settings;
//        using Base_T::edge_count;
//        using Base_T::Settings;
//        using Base_T::random_engine;
//        using Base_T::normal_dist;
//
//        Sampler() : Base_T() {}
//
//        virtual ~Sampler()
//        {
//            safe_free(x_buffer);
//            safe_free(y_buffer);
//            safe_free(p_buffer);
//            safe_free(r);
//            safe_free(rho);
//        }
//
//        explicit Sampler(
//            const Int edge_count_,
//            const Setting_T settings_ = Setting_T()
//        )
//        :   Base_T( edge_count_, settings_ )
//        {
//            safe_alloc ( x_buffer,  edge_count    * AmbDim );
//            safe_alloc ( y_buffer,  edge_count    * AmbDim );
//            safe_alloc ( p_buffer, (edge_count+1) * AmbDim );
//
//            safe_alloc (r,   edge_count );
//            fill_buffer(r,   edge_count, one / edge_count);
//
//            safe_alloc (rho, edge_count );
//            fill_buffer(rho, edge_count, one / edge_count);
//
//
//            total_r_inv = one;
//        }
//
//        explicit Sampler(
//            const Real * restrict const r_in,
//            const Real * restrict const rho_in,
//            const Int edge_count_,
//            const Setting_T settings_ = Setting_T()
//        )
//        :   Base_T( edge_count_, settings_ )
//        {
//            safe_alloc( x_buffer,  edge_count    * AmbDim );
//            safe_alloc( y_buffer,  edge_count    * AmbDim );
//            safe_alloc( p_buffer, (edge_count+1) * AmbDim );
//
//            safe_alloc( r,   edge_count );
//            safe_alloc( rho, edge_count );
//
//            ReadEdgeLengths(r_in);
//            ReadRho(rho_in);
//        }
//
//
//
////        // Copy constructor
////        Sampler( const Sampler & other )
////        :   Base_T( other.edge_count, other.settings )
////        ,   x(other.x)
////        ,   y(other.y)
////        ,   p(other.p)
////        ,   r(other.r)
////        ,   rho(other.rho)
////        ,   total_r_inv(other.total_r_inv)
////        ,   w(other.w)
////        ,   F(other.F)
////        ,   DF(other.DF)
////        ,   L(other.L)
////        ,   u(other.u)
////        ,   iter(other.iter)
////        ,   residual(other.residual)
////        {}
//
////
////        friend void swap(Sampler &A, Sampler &B) noexcept
////        {
////            // see https://stackoverflow.com/questions/5695548/public-friend-swap-member-function for details
////            using std::swap;
////
////            swap(A.edge_count,B.edge_count);
////            swap(A.x,B.x);
////            swap(A.y,B.y);
////            swap(A.p,B.p);
////            swap(A.r,B.r);
////            swap(A.rho,B.rho);
////            swap(A.total_r_inv,B.total_r_inv);
////
////            std::swap_ranges(&A.w[0],    &A.w[AmbDim],            &B.w[0]   );
////            std::swap_ranges(&A.z[0],    &A.z[AmbDim],            &B.z[0]   );
////            std::swap_ranges(&A.F[0],    &A.F[AmbDim],            &B.F[0]   );
////            std::swap_ranges(&A.A(0,0),  &A.A(AmbDim-1,AmbDim),   &B.A(0,0) );
////            std::swap_ranges(&A.u[0],    &A.u[AmbDim],            &B.u[0]   );
////
////            swap(A.settings,       B.settings     );
////            swap(A.iter,           B.iter          );
////            swap(A.residual,       B.residual      );
////        }
////
////        // Copy assignment operator
////        Sampler & operator=(Sampler other)
////        {
////            // copy-and-swap idiom
////            // see https://stackoverflow.com/a/3279550/8248900 for details
////            swap(*this, other);
////
////            return *this;
////        }
////
////        /* Move constructor */
////        Sampler( Sampler && other ) noexcept
////        :   Sampler()
////        {
////            swap(*this, other);
////        }
//
//
//    protected:
//
//        mutable Real restrict * x_buffer = nullptr;
//        mutable Real restrict * y_buffer = nullptr;
//        mutable Real restrict * p_buffer = nullptr;
//
//        Real restrict * r   = nullptr;
//        Real restrict * rho = nullptr;
//
//        Real total_r_inv = one;
//
//        Real DF [AmbDim][AmbDim];  // nabla F(0) with respect to measure ys
//        Real L  [AmbDim][AmbDim];  // storing Cholesky factor.
//
//        Real w [AmbDim];           // current point in hyperbolic space.
//        Real F [AmbDim];           // right hand side of Newton iteration.
//        Real u [AmbDim];           // update direction
//        Real z [AmbDim];           // Multiple purpose buffer.
//
//        ShiftMap<AmbDim,Real,Int> S;
//
//        Int iter = 0;
//
//        Real squared_residual = 1;
//        Real         residual = 1;
//
//
//        Real edge_space_sampling_weight = 0;
//        Real edge_quotient_space_sampling_correction = 0;
//
//
//        Real lambda_min = eps;
//        Real q = one;
//        Real errorestimator = infty;
//
//        bool linesearchQ = true;    // Toggle line search.
//        bool succeededQ  = false;   // Whether algorithm has succeded.
//        bool continueQ   = true;    // Whether to continue with the main loop.
//        bool ArmijoQ     = false;   // Whether Armijo condition was met last time we checked.
//
//        static constexpr Real zero              = 0;
//        static constexpr Real half              = 0.5;
//        static constexpr Real one               = 1;
//        static constexpr Real two               = 2;
//        static constexpr Real three             = 3;
//        static constexpr Real four              = 4;
//        static constexpr Real eps               = std::numeric_limits<Real>::min();
//        static constexpr Real infty             = std::numeric_limits<Real>::max();
//        static constexpr Real small_one         = 1 - 16 * eps;
//        static constexpr Real big_one           = 1 + 16 * eps;
//        static constexpr Real g_factor          = 4;
//        static constexpr Real g_factor_inv      = one/g_factor;
//        static constexpr Real norm_threshold    = 0.99 * 0.99 + 16 * eps;
//        static constexpr Real two_pi            = static_cast<Real>(2 * M_PI);
//
//    public:
//
//        virtual void Optimize() override
//        {
//            const Int max_iter = settings.max_iter;
//
//            iter = 0;
//
//            Shift();
//
//            DifferentialAndHessian_Hyperbolic();
//
//            SearchDirection_Hyperbolic();
//
//            // The main loop.
//            while( ( iter < max_iter ) && continueQ )
//            {
//                ++iter;
//
//                LineSearch_Hyperbolic_Potential();
//
//                DifferentialAndHessian_Hyperbolic();
//
//                SearchDirection_Hyperbolic();
//            }
//        }
//
//
//        void ComputeEdgeSpaceSamplingWeight()
//        {
//            edge_space_sampling_weight = S.EdgeSpaceSamplingWeight(
//                x_buffer, &w[0], y_buffer, r, rho, edge_count
//            );
//        }
//
//        virtual Real EdgeSpaceSamplingWeight() const override
//        {
//            return edge_space_sampling_weight;
//        }
//
//        void ComputeEdgeQuotientSpaceSamplingCorrection()
//        {
//            if( AmbDim == 2)
//            {
//                edge_quotient_space_sampling_correction = one;
//                return;
//            }
//
//            Eigen::Matrix<Real,AmbDim,AmbDim> Sigma;
//
//            // We fill only the lower triangle of Sigma, because that's the only thing that Eigen' selfadjoint eigensolver needs.
//            // Recall that Eigen matrices are column-major by default.
//
//            {
//                const Real rho_squared = rho[0] * rho[0];
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    const Real factor = rho_squared * y(0,i);
//
//                    for( Int j = i; j < AmbDim; ++j )
//                    {
//                        Sigma(j,i) = factor * y(0,j);
//                    }
//                }
//            }
//
//            for( Int k = 0; k < edge_count; ++k )
//            {
//                const Real rho_squared = rho[k] * rho[k];
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    const Real factor = rho_squared * y(k,i);
//
//                    for( Int j = i; j < AmbDim; ++j )
//                    {
//                        Sigma(j,i) += factor * y(k,j);
//                    }
//                }
//            }
//
//            // Eigen needs only the lower triangular part. So need not symmetrize.
//
////            for( Int i = 0; i < AmbDim; ++i )
////            {
////                for( Int j = 0; j < i; ++j )
////                {
////                    Sigma(j,i) = Sigma(i,j);
////                }
////            }
//
//            if( AmbDim == 3)
//            {
//                // Exploiting that
//                //      (lambda[0] + lambda[1]) * (lambda[0] + lambda[2]) * (lambda[1] + lambda[2])
//                //      =
//                //      ( tr(Sigma*Sigma) - tr(Sigma)*tr(Sigma) ) *  tr(Sigma)/2 - det(Sigma)
//                //  Thus, it can be expressed by as third-order polynomial in the entries of the matrix.
//
//                const Real S_00 = Sigma(0,0)*Sigma(0,0);
//                const Real S_11 = Sigma(1,1)*Sigma(1,1);
//                const Real S_22 = Sigma(2,2)*Sigma(2,2);
//
//                const Real S_10 = Sigma(1,0)*Sigma(1,0);
//                const Real S_20 = Sigma(2,0)*Sigma(2,0);
//                const Real S_21 = Sigma(2,1)*Sigma(2,1);
//
//                const Real det = std::abs(
//                      Sigma(0,0) * ( S_11 + S_22 - S_10 - S_20 )
//                    + Sigma(1,1) * ( S_00 + S_22 - S_10 - S_21 )
//                    + Sigma(2,2) * ( S_00 + S_11 - S_20 - S_21 )
//                    + two * (Sigma(0,0)*Sigma(1,1)*Sigma(2,2) - Sigma(1,0)*Sigma(2,0)*Sigma(2,1))
//                );
//                edge_quotient_space_sampling_correction = one / std::sqrt(det);
//                return;
//            }
//
//
//            Eigen::SelfAdjointEigenSolver< Eigen::Matrix<Real,AmbDim,AmbDim> > eigs;
//
//            eigs.compute(Sigma);
//
//            auto & lambda = eigs.eigenvalues();
//
//            Real det = one;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                for( Int j = i+1; j < AmbDim; ++j )
//                {
//                    det *= (lambda(i)+lambda(j));
//                }
//            }
//
//            edge_quotient_space_sampling_correction = one / std::sqrt(det);
//        }
//
//        virtual Real EdgeQuotientSpaceSamplingCorrection() const override
//        {
//            return edge_quotient_space_sampling_correction;
//        }
//
//        virtual Real EdgeQuotientSpaceSamplingWeight() const override
//        {
//            return EdgeSpaceSamplingWeight() * EdgeQuotientSpaceSamplingCorrection();
//        }
//
//
//    protected:
//
//
//
////        Real & x( const Int k, const Int i )
////        {
////            return x_buffer[AmbDim*k+i];
////        }
////
////        const Real & x( const Int k, const Int i ) const
////        {
////            return x_buffer[AmbDim*k+i];
////        }
////
////        Real & y( const Int k, const Int i )
////        {
////            return y_buffer[AmbDim*k+i];
////        }
////
////        const Real & y( const Int k, const Int i ) const
////        {
////            return y_buffer[AmbDim*k+i];
////        }
////
////        Real & p( const Int k, const Int i )
////        {
////            return p_buffer[AmbDim*k+i];
////        }
////
////        const Real & p( const Int k, const Int i ) const
////        {
////            return p_buffer[AmbDim*k+i];
////        }
//
//        force_inline Real Dot_ww() const
//        {
//            Real ww = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                ww += w[i] * w[i];
//            }
//            return ww;
//        }
//
//        force_inline Real Dot_zz() const
//        {
//            Real zz = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                zz += z[i] * z[i];
//            }
//            return zz;
//        }
//
//        force_inline Real Dot_wz() const
//        {
//            Real wz = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                wz += w[i] * z[i];
//            }
//            return wz;
//        }
//
//        force_inline Real Dot_uu() const
//        {
//            Real uu = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                uu += u[i] * u[i];
//            }
//            return uu;
//        }
//
//        force_inline Real Dot_Fu() const
//        {
//            Real Fu = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                Fu += F[i] * u[i];
//            }
//            return Fu;
//        }
//
//        force_inline Real Dot_FF() const
//        {
//            Real FF = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                FF += F[i] * F[i];
//            }
//            return FF;
//        }
//
//        force_inline Real Dot_FDFu() const
//        {
//            Real result = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                Real z_i = 0;
//                for( Int j = 0; j < i; ++j )
//                {
//                    z_i += DF[j][i] * u[j];
//                }
//                for( Int j = i; j < AmbDim; ++j )
//                {
//                    z_i += DF[i][j] * u[j];
//                }
//
//                result += F[i] * z_i;
//
//            }
//
//            return result;
//        }
//
//        void L_Cholesky()
//        {
//            for( Int k = 0; k < AmbDim; ++k )
//            {
//                const Real a = L[k][k] = std::sqrt(L[k][k]);
//                const Real ainv = one/a;
//
//                for( Int j = k+1; j < AmbDim; ++j )
//                {
//                    L[k][j] *= ainv;
//                }
//
//                for( Int i = k+1; i < AmbDim; ++i )
//                {
//                    for( Int j = i; j < AmbDim; ++j )
//                    {
//                        L[i][j] -= L[k][i] * L[k][j];
//                    }
//                }
//            }
//        }
//
//        void L_CholeskySolve_Fu()
//        {
//            copy_buffer( &F[0], &u[0], AmbDim);
//
//            // Now solve in-place.
//
//            // Lower triangular back substitution
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                for( Int j = 0; j < i; ++j )
//                {
//                    u[i] -= L[j][i] * u[j];
//                }
//                u[i] /= L[i][i];
//            }
//
//            // Upper triangular back substitution
//            for( Int i = AmbDim; i --> 0; )
//            {
//                for( Int j = i+1; j < AmbDim; ++j )
//                {
//                    u[i] -= L[i][j] * u[j];
//                }
//                u[i] /= L[i][i];
//            }
//        }
//
//        Real DF_SmallestEigenValue() const
//        {
//            if( AmbDim == 2)
//            {
//                Real lambda_min = half * (
//                    DF[0][0] + DF[1][1]
//                    - std::sqrt(
//                        std::abs(
//                            (DF[0][0]-DF[1][1])*(DF[0][0]-DF[1][1]) + four * DF[0][1]*DF[0][1]
//                        )
//                    )
//                );
//
//                return lambda_min;
//            }
//
//            if( AmbDim == 3)
//            {
//                Real lambda_min;
//
//                const Real p1 = DF[0][1]*DF[0][1] + DF[0][2]*DF[0][2] + DF[1][2]*DF[1][2];
//
//                if( std::sqrt(p1) < eps * std::sqrt( DF[0][0]*DF[0][0] + DF[1][1]*DF[1][1] + DF[2][2]*DF[2][2]) )
//                {
//                    // DF is diagonal
//                    lambda_min = std::min( DF[0][0], std::min(DF[1][1],DF[2][2]) );
//                }
//                else
//                {
//                    const Real q         = ( DF[0][0] + DF[1][1] + DF[2][2] ) / three;
//                    const Real delta [3] = { DF[0][0]-q, DF[1][1]-q, DF[2][2]-q } ;
//                    const Real p2   = delta[0]*delta[0] + delta[1]*delta[1] + delta[2]*delta[2] + two*p1;
//                    const Real p    = std::sqrt( p2 / static_cast<Real>(6) );
//                    const Real pinv = one/p;
//                    const Real b11  = delta[0] * pinv;
//                    const Real b22  = delta[1] * pinv;
//                    const Real b33  = delta[2] * pinv;
//                    const Real b12  = DF[0][1] * pinv;
//                    const Real b13  = DF[0][2] * pinv;
//                    const Real b23  = DF[1][2] * pinv;
//
//                    const Real r = half * (two * b12 * b23 * b13 - b11 * b23 * b23 - b12 *b12 * b33 + b11 * b22 * b33 - b13 *b13 * b22);
//
//
//                    const Real phi = ( r <= -one )
//                        ? ( static_cast<Real>(M_PI) / three )
//                        : ( ( r >= one ) ? zero : acos(r) / three );
//
//                    // The eigenvalues are ordered this way: eig2 <= eig1 <= eig0.
//
////                    Real eig0 = q + two * p * cos( phi );
////                    Real eig2 = q + two * p * cos( phi + two * M_PI/ three );
////                    Real eig1 = three * q - eig0 - eig2;
//
//                    lambda_min = q + two * p * cos( phi + two * M_PI/ three );
//                }
//
//                return lambda_min;
//            }
//
//            using Matrix_T = Eigen::Matrix<Real,AmbDim,AmbDim>;
//
//            Matrix_T Sigma (&DF[0][0]);
//
//            Eigen::SelfAdjointEigenSolver<Matrix_T> eigs;
//
//            eigs.compute(Sigma);
//
//            return eigs.eigenvalues()[0];
//        }
//
//        Real Potential()
//        {
//            Real value = 0;
//
//            const Real zz = Dot_zz();
//
//            const Real a = big_one + zz;
//            const Real c = (big_one-zz);
//
//            const Real b = one/c;
//
//            for( Int k = 0; k < edge_count; ++k )
//            {
//                Real yz2 = 0;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    yz2 += y(k,i) * z[i];
//                }
//
//                value += r[k] * std::log(std::abs( (a - two * yz2) * b ) );
//            }
//
//            return value * total_r_inv;
//        }
//
//
//        void LineSearch_Hyperbolic_Residual()
//        {
//            // 2 F(0)^T.DF(0).u is the derivative of w\mapsto F(w)^T.F(w) at w = 0.
//
//            const Real slope = two * Dot_FDFu();
//
//            Real tau = one;
//
//            const Real u_norm = std::sqrt(Dot_uu());
//
//            // exponential map shooting from 0 to tau * u.
//            {
//                const Real scale = tau * tanhc(tau * u_norm);
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    z[i] = scale * u[i];
//                }
//            }
//
//            // Shift the point z along -w to get new updated point w .
//            InverseShift();
//
//            // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
//            Shift();
//
//            const Real squared_residual_at_0 = squared_residual;
//
//            DifferentialAndHessian_Hyperbolic();
//
//            if( linesearchQ )
//            {
//
//                Int backtrackings = 0;
//
//                // Armijo condition for squared residual.
//
//                ArmijoQ = squared_residual - squared_residual_at_0 - settings.Armijo_slope_factor * tau * slope < static_cast<Real>(0);
//
//                while( !ArmijoQ && (backtrackings < settings.max_backtrackings) )
//                {
//                    ++backtrackings;
//
//                    // Estimate step size from quadratic fit if applicable.
//
//                    const Real tau_1 = settings.Armijo_shrink_factor * tau;
//                    const Real tau_2 = - half * settings.Armijo_slope_factor * tau * tau * slope / ( squared_residual  - squared_residual_at_0 - tau * slope );
//
//                    tau = std::max( tau_1, tau_2 );
//
//                    {
//                        const Real scale = tau * tanhc(tau * u_norm);
//
//                        for( Int i = 0; i < AmbDim; ++i )
//                        {
//                            z[i] = scale * u[i];
//                        }
//                    }
//
//                    // Shift the point z along -w to get new updated point w .
//                    InverseShift();
//
//                    // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
//                    Shift();
//
//                    DifferentialAndHessian_Hyperbolic();
//
//
//                    ArmijoQ = squared_residual - squared_residual_at_0 - settings.Armijo_slope_factor * tau * slope < 0;
//
//                }
//            }
//        }
//
//        void LineSearch_Hyperbolic_Potential()
//        {
//            Real tau = one;
//
//            const Real u_norm = std::sqrt(Dot_uu());
//
//            // exponential map shooting from 0 to tau * u.
//            {
//                const Real scale = tau * tanhc(tau * u_norm);
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    z[i] = scale * u[i];
//                }
//            }
//
//            if( linesearchQ )
//            {
//                //Linesearch with potential as merit function.
//
//                const Real gamma = settings.Armijo_shrink_factor;
//
//                const Real sigma = settings.Armijo_slope_factor;
//
//                const Real Dphi_0 = g_factor * Dot_Fu();
//
//                Int backtrackings = 0;
//
//                // Compute potential and check Armijo condition.
//
////                const Real phi_0 = Potential(o);
//
//                Real phi_tau = Potential();
//
//                ArmijoQ = phi_tau /*- phi_0*/ - sigma * tau * Dphi_0 < 0; // phi_0 == 0
//
//
//                while( !ArmijoQ && (backtrackings < settings.max_backtrackings) )
//                {
//                    ++backtrackings;
//
//                    const Real tau_1 = gamma * tau;
//
//                    // Estimate step size from quadratic fit if applicable.
//                    const Real tau_2 = - half * sigma * tau * tau * Dphi_0 / ( phi_tau /*- phi_0*/ - tau * Dphi_0 );
//
//                    tau = std::max( tau_1, tau_2 );
//
//                    {
//                        const Real scale = tau * tanhc(tau * u_norm);
//
//                        for( Int i = 0; i < AmbDim; ++i )
//                        {
//                            z[i] = scale * u[i];
//                        }
//                    }
//
//                    phi_tau = Potential();
//
//                    ArmijoQ = phi_tau  /*- phi_0*/ - sigma * tau * Dphi_0 < 0; // phi_0 == 0
//                }
//            }
//
//            // Shift the point z along -w to get new updated point w .
//            InverseShift();
//
//            // Shift the input measure along w to 0 to simplify gradient, Hessian, and update computation .
//            Shift();
//        }
//
//        void DifferentialAndHessian_Hyperbolic()
//        {
//            // CAUTION: We use a different sign convention as in the paper!
//            // Assemble  F = -1/2 y * r.
//            // Assemble DF = nabla F + regulatization:
//            // DF_{ij} = \delta_{ij} - \sum_k x_{k,i} x_{k,j} \r_k.
//
//            // Overwrite F and DF by contribution of first element of y.
//            {
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    const Real factor = r[0] * y(0,i);
//
//                    F[i] = - factor;
//
//                    for( Int j = i; j < AmbDim; ++j )
//                    {
//                        DF[i][j] = - factor * y(0,j);
//                    }
//                }
//            }
//
//            // Add-in the contributions of the other elements of y.
//            for( Int k = 1; k < edge_count; ++k )
//            {
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    const Real factor = r[k] * y(k,i);
//
//                    F[i] -= factor;
//
//                    for( Int j = i; j < AmbDim; ++j )
//                    {
//                        DF[i][j] -= factor * y(k,j);
//                    }
//                }
//            }
//
//            // Normalize for case that the weights in r do not sum to 1.
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                F[i] *= total_r_inv;
//
//                for( Int j = i; j < AmbDim; ++j )
//                {
//                    DF[i][j] *= total_r_inv;
//                }
//            }
//
//            squared_residual = Dot_FF();
//
//            residual = std::sqrt( squared_residual );
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                F[i] *= half;
//                // We better add the identity for precision reasons here -- after the other computations.
//                DF[i][i] += one;
//            }
//        }
//
//        void SearchDirection_Hyperbolic()
//        {
//            // Make decisions whether to continue.
//            if( residual < static_cast<Real>(100.) * settings.tolerance )
//            {
//                // We have to compute eigenvalue _before_ we add the regularization.
//
//                lambda_min = DF_SmallestEigenValue();
//
//                q = four * residual / (lambda_min * lambda_min);
//
//                if( q < one )
//                {
//                    //Kantorovich condition satisfied; this allows to compute an error estimator.
//                    errorestimator = half * lambda_min * q;
//                    //And we should deactivate line search. Otherwise, we may run into precision issues.
//                    linesearchQ = false;
//                    continueQ = (errorestimator >settings.tolerance);
//                    succeededQ = !continueQ;
//                }
//                else
//                {
//                    errorestimator = infty;
//                    linesearchQ = settings.Armijo_slope_factor > zero;
//                    //There is no way to reduce the residual below machine epsilon. If the algorithm reaches here, the problem is probably too ill-conditioned to be solved in machine precision.
//                    continueQ = residual > settings.give_up_tolerance;
//                }
//            }
//            else
//            {
//                q = big_one;
//                lambda_min = eps;
//                errorestimator = infty;
//                linesearchQ = settings.Armijo_slope_factor > zero;
//                continueQ = residual>std::max( settings.give_up_tolerance, settings.tolerance );
//            }
//
//            const Real c = settings.regularization * squared_residual;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                for( Int j = i; j < AmbDim; ++j )
//                {
//                    L[i][j] = DF[i][j] + static_cast<Real>(i==j) * c;
//                }
//            }
//
//            L_Cholesky();
//
//            L_CholeskySolve_Fu();
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                u[i] *= -one;
//            }
//        }
//
//        void Gradient_Hyperbolic()
//        {
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                u[i] = -g_factor_inv * F[i];
//            }
//        }
//
//        void Gradient_Planar()
//        {
//            // The factor 2 is here to reproduce the Abikoff-Ye algorithm (in the absence of linesearch.)
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                u[i] = -two * F[i];
//            }
//        }
//
//    public:
//
//        void InverseShift()
//        {
//            // Shift point w along -z.
//
//            const Real wz2 = two * Dot_wz();
//            const Real ww  = Dot_ww();
//            const Real zz  = Dot_zz();
//
//            const Real a = one - ww;
//            const Real b = one + zz + wz2;
//            const Real c = big_one + wz2 + ww * zz;
//            const Real d = one / c;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                w[i] = ( a * z[i] + b * w[i] ) * d;
//            }
//        }
//
//        void Shift()
//        {
//            S.Shift( x_buffer, &w[0], y_buffer, edge_count, one );
//        }
//
//    public:
//
////
////        virtual const SpherePoints_T & InitialEdgeCoordinates() const override
////        {
////            return x;
////        }
//
//        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, bool normalize = true ) override
//        {
//            copy_buffer(x_in, x_buffer, edge_count * AmbDim );
//
//            if( normalize )
//            {
//                for( Int k = 0; k < edge_count; ++k )
//                {
//                    Real r2 = 0;
//
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        r2 += x(k,i) * x(k,i);
//                    }
//
//                    const Real scale = one/std::sqrt(r2);
//
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        x(k,i) *= scale;
//                    }
//                }
//            }
//        }
//
//        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, const Int k, bool normalize = true ) override
//        {
//            ReadInitialEdgeCoordinates( &x_in[ AmbDim * edge_count * k], normalize);
//        }
//
//        virtual void WriteInitialEdgeCoordinates( Real * x_out ) const override
//        {
//            copy_buffer( x_buffer, x_out, edge_count * AmbDim );
//        }
//
//        virtual void WriteInitialEdgeCoordinates( Real * x_out, const Int k ) const override
//        {
//            WriteInitialEdgeCoordinates( &x_out[ AmbDim * edge_count * k ]);
//        }
//
//        virtual const Real * EdgeCoordinates() const override
//        {
//            return y_buffer;
//        }
//
//        virtual void ReadEdgeCoordinates( const Real * const y_in ) override
//        {
//            copy_buffer( y_in, y_buffer, edge_count * AmbDim );
//        }
//
//        virtual void ReadEdgeCoordinates( const Real * const y_in, const Int k ) override
//        {
//            ReadEdgeCoordinates( &y_in[ AmbDim * edge_count * k ]);
//        }
//
//        virtual void WriteEdgeCoordinates( Real * y_out ) const override
//        {
//            copy_buffer( y_buffer, y_out, edge_count * AmbDim );
//        }
//
//        virtual void WriteEdgeCoordinates( Real * y_out, const Int k ) const override
//        {
//            WriteEdgeCoordinates( &y_out[ AmbDim * edge_count * k ]);
//        }
//
//
//
//        virtual const Real * SpaceCoordinates() const override
//        {
//            return p_buffer;
//        }
//
//        virtual void WriteSpaceCoordinates( Real * p_out ) const  override
//        {
//            copy_buffer( p_buffer, p_out, (edge_count+1)*AmbDim );
//        }
//
//        virtual void WriteSpaceCoordinates( Real * p_out, const Int k ) const override
//        {
//            WriteSpaceCoordinates( &p_out[ (edge_count+1) * AmbDim * k ]);
//        }
//
//        void ComputeSpaceCoordinates() override
//        {
//            //Caution: This gives only have the weight to the end vertices of the chain.
//            //Thus this is only really the barycenter, if the chain is closed!
//
//            Real barycenter        [AmbDim] = {};
//            Real point_accumulator [AmbDim] = {};
//
//            for( Int k = 0; k < edge_count; ++k )
//            {
//                const Real r_k = r[k];
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    const Real offset = r_k * y(k,i);
//
//                    barycenter[i] += (point_accumulator[i] + half * offset);
//
//                    point_accumulator[i] += offset;
//                }
//            }
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                p(0,i) = -barycenter[i]/edge_count;
//            }
//
//            for( Int k = 0; k < edge_count; ++k )
//            {
//                const Real r_k = r[k];
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    p(k+1,i) = p(k,i) + r_k * y(k,i);
//                }
//            }
//        }
//
//        virtual const Real * EdgeLengths() const override
//        {
//            return r;
//        }
//
//        virtual void ReadEdgeLengths( const Real * const r_in ) override
//        {
//            copy_buffer(r_in, r, edge_count);
//
//            Real sum = 0;
//
//            for( Int k = 0; k < edge_count; ++k )
//            {
//                sum += r[k];
//            }
//
//            total_r_inv = one / sum;
//        }
//
//        virtual const Real * Rho() const override
//        {
//            return rho;
//        }
//
//        virtual void ReadRho( const Real * const rho_in ) override
//        {
//            copy_buffer(rho_in, rho, edge_count);
//        }
//
//
//        virtual void ComputeShiftVector() override
//        {
//            Real w_ [AmbDim] = {};
//
//            for( Int k = 0; k < edge_count; ++k )
//            {
//                const Real r_k = r[k];
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    w_[i] += x(k,i) * r_k;
//                }
//            }
//
//            // Normalize in that case that r does not sum up to 1.
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                w[i] = w_[i] * total_r_inv;
//            }
//        }
//
//        virtual void ReadShiftVector( const Real * const w_in ) override
//        {
//            copy_buffer(&w_in[0], &w[0], AmbDim);
//
//            // Use Euclidean barycenter as initial guess if the supplied initial guess does not make sense.
//            if( Dot_ww() > small_one )
//            {
//                ComputeShiftVector();
//            }
//        }
//
//        virtual void ReadShiftVector( const Real * const w_in, const Int k ) override
//        {
//            ReadShiftVector( &w_in[ AmbDim * k ] );
//        }
//
//        virtual void WriteShiftVector( Real * w_out ) const override
//        {
//            copy_buffer(&w[0],&w_out[0],AmbDim);
//        }
//
//        virtual void WriteShiftVector( Real * w_out, const Int k ) const override
//        {
//            copy_buffer(&w[0],&w_out[ AmbDim * k],AmbDim);
//        }
//
//        virtual Real Residual() const override
//        {
//            return residual;
//        }
//
//        virtual Real ErrorEstimator() const override
//        {
//            return errorestimator;
//        }
//
//        virtual Int IterationCount() const override
//        {
//            return iter;
//        }
//
//        virtual Int MaxIterationCount() const override
//        {
//            return settings.max_iter;
//        }
//
//
//    public:
//
//        virtual void OptimizeBatch(
//            const Real * const x_in,
//                  Real *       w_out,
//                  Real *       y_out,
//            const Int sample_count,
//            const Int thread_count = 1,
//            bool normalize = true
//        ) override
//        {
//            ptic(ClassName()+"OptimizeBatch");
//
//            JobPointers<Int> job_ptr ( sample_count, thread_count );
//
//            #pragma omp parallel for num_threads( thread_count )
//            for( Int thread = 0; thread < thread_count; ++thread )
//            {
//                const Int k_begin = job_ptr[thread];
//                const Int k_end   = job_ptr[thread+1];
//
//                Sampler W( edge_count, settings );
//
//                W.ReadEdgeLengths( r );
//
//                for( Int k = k_begin; k < k_end; ++k )
//                {
//                    W.ReadInitialEdgeCoordinates( x_in, k, normalize );
//
//                    W.ComputeShiftVector();
//
//                    W.Optimize();
//
//                    W.WriteShiftVector( w_out, k );
//
//                    W.WriteEdgeCoordinates( y_out, k );
//                }
//            }
//
//            ptoc(ClassName()+"OptimizeBatch");
//        }
//
//        virtual void RandomClosedPolygons(
//                  Real * const restrict x_out,
//                  Real * const restrict w_out,
//                  Real * const restrict y_out,
//                  Real * const restrict K_edge_space,
//                  Real * const restrict K_edge_quotient_space,
//            const Int sample_count,
//            const Int thread_count = 1
//        ) const override
//        {
//            ptic(ClassName()+"RandomClosedPolygons");
//
//            JobPointers<Int> job_ptr ( sample_count, thread_count );
//
////            valprint( "dimension   ", AmbDim       );
////            valprint( "edge_count  ", edge_count   );
////            valprint( "sample_count", sample_count );
////            valprint( "thread_count", thread_count );
//
//            #pragma omp parallel for num_threads( thread_count )
//            for( Int thread = 0; thread < thread_count; ++thread )
//            {
//                const Int k_begin = job_ptr[thread  ];
//                const Int k_end   = job_ptr[thread+1];
//
//                Sampler W( edge_count, settings );
//
//                W.ReadEdgeLengths( r );
//                W.ReadRho( rho );
//
//                for( Int k = k_begin; k < k_end; ++k )
//                {
//                    W.RandomizeInitialEdgeCoordinates();
//
//                    W.WriteInitialEdgeCoordinates(x_out, k);
//
//                    W.ComputeShiftVector();
//
//                    W.Optimize();
//
//                    W.WriteShiftVector(w_out, k);
//
//                    W.WriteEdgeCoordinates(y_out, k);
//
//                    W.ComputeEdgeSpaceSamplingWeight();
//
//                    W.ComputeEdgeQuotientSpaceSamplingCorrection();
//
//                    K_edge_space[k] = W.EdgeSpaceSamplingWeight();
//
//                    K_edge_quotient_space[k] = W.EdgeQuotientSpaceSamplingWeight();
//                }
//            }
//
//            ptoc(ClassName()+"::RandomClosedPolygons");
//        }
//
//
//
//
//
//
//        // moments: A 3D-array of size 3 x fun_count x bin_count. Entry moments(i,j,k) will store the sampled weighted k-th moment of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i (see above).
//        // ranges: Specify the range for binning: For j-th function in F_list, the range from ranges(j,0) to ranges(j,1) will be devided into bin_count bins. The user is supposed to provide meaningful ranges. Some rough guess might be obtained by calling the random variables on the prepared Sampler_T C.
//        void Sample_Binned(
//            Real * restrict bins_out,
//            const Int bin_count_,
//            Real * restrict moments_out,
//            const Int moment_count_,
//            const Real * restrict ranges,
//            const std::vector< std::unique_ptr<RandomVariable_T> > & F_list_,
//            const Int sample_count,
//            const Int thread_count = 1
//        ) const
//        {
//            ptic(ClassName()+"Sample_Binned (polymorphic)");
//
//            const Int fun_count = static_cast<Int>(F_list_.size());
//
//            const Int moment_count = std::max( static_cast<Int>(3), moment_count_ );
//
//            const Int bin_count = std::max( bin_count_, static_cast<Int>(1) );
//
//            JobPointers<Int> job_ptr ( sample_count, thread_count );
//
//            valprint( "dimension   ", AmbDim       );
//            valprint( "edge_count  ", edge_count );
//            valprint( "sample_count", sample_count );
//            valprint( "fun_count   ", fun_count );
//            valprint( "bin_count   ", bin_count );
//            valprint( "moment_count", moment_count );
//            valprint( "thread_count", thread_count );
//
//
//            Tensor3<Real,Int> bins_global   ( bins_out,    3, fun_count, bin_count    );
//            Tensor3<Real,Int> moments_global( moments_out, 3, fun_count, moment_count );
//            Tensor1<Real,Int> factor        (                 fun_count               );
//
//            print("Sampling the following random variables:");
//            for( Int i = 0; i < fun_count; ++ i )
//            {
//                const size_t i_ = static_cast<size_t>(i);
//                factor(i) = static_cast<Real>(bin_count) / ( ranges[2*i+1] - ranges[2*i+0] );
//
//                print("    " + F_list_[i_]->Tag());
//            }
//
//            const Int lower = static_cast<Int>(0);
//            const Int upper = static_cast<Int>(bin_count-1);
//
//            #pragma omp parallel for num_threads( thread_count )
//            for( Int thread = 0; thread < thread_count; ++thread )
//            {
//                const Int repetitions = (job_ptr[thread+1] - job_ptr[thread]);
//
//                Sampler W( edge_count, settings );
//
//                W.ReadEdgeLengths( r );
//                W.ReadRho( rho );
//
//                std::vector< std::unique_ptr<RandomVariable_T> > F_list (fun_count);
//
//                for( Int i = 0; i < fun_count; ++ i )
//                {
//                    F_list[i] = F_list_[static_cast<size_t>(i)]->Clone();
//                }
//
//                Tensor3<Real,Int> bins_local   ( 3, fun_count, bin_count,    zero );
//                Tensor3<Real,Int> moments_local( 3, fun_count, moment_count, zero );
//
//                for( Int k = 0; k < repetitions; ++k )
//                {
//                    W.RandomizeInitialEdgeCoordinates();
//
//                    W.ComputeShiftVector();
//
//                    W.Optimize();
//
//                    W.ComputeSpaceCoordinates();
//
//                    W.ComputeEdgeSpaceSamplingWeight();
//
//                    W.ComputeEdgeQuotientSpaceSamplingCorrection();
//
//                    const Real K = W.EdgeSpaceSamplingWeight();
//
//                    const Real K_quot = W.EdgeQuotientSpaceSamplingWeight();
//
//                    for( Int i = 0; i < fun_count; ++i )
//                    {
//                        auto & F = *F_list[static_cast<size_t>(i)];
//
//                        const Real val = F(W);
//
//                        Real values [3] = {static_cast<Real>(1),K,K_quot};
//
////                        const Int bin_idx = std::clamp(
////                           static_cast<Int>(std::floor( factor[i] * (val - ranges(i,0)) )),
////                           lower,
////                           upper
////                        );
//
//                        const Int bin_idx = static_cast<Int>(std::floor( factor[i] * (val - ranges[2*i]) ));
//
//                        if( (bin_idx <= upper) && (bin_idx >= lower) )
//                        {
//                            bins_local(0,i,bin_idx) += static_cast<Real>(1);
//                            bins_local(1,i,bin_idx) += K;
//                            bins_local(2,i,bin_idx) += K_quot;
//                        }
//
//                        moments_local(0,i,0) += values[0];
//                        moments_local(1,i,0) += values[1];
//                        moments_local(2,i,0) += values[2];
//
//                        for( Int j = 1; j < moment_count; ++j )
//                        {
//                            values[0] *= val;
//                            values[1] *= val;
//                            values[2] *= val;
//                            moments_local(0,i,j) += values[0];
//                            moments_local(1,i,j) += values[1];
//                            moments_local(2,i,j) += values[2];
//                        }
//                    }
//                }
//
//                #pragma omp critical
//                {
////                    valprint("thread",thread);
////                    print( "data_local = " + data_local.ToString() );
//                    for( Int l = 0; l < 3; ++l )
//                    {
//                        for( Int i = 0; i < fun_count; ++i )
//                        {
//                            for( Int j = 0; j < bin_count; ++j )
//                            {
//                                bins_global(l,i,j) += bins_local(l,i,j);
//                            }
//
//                            for( Int j = 0; j < moment_count; ++j )
//                            {
//                                moments_global(l,i,j) += moments_local(l,i,j);
//                            }
//                        }
//                    }
//                }
//            }
//
//            bins_global.Write( bins_out );
//            moments_global.Write( moments_out );
//
//            ptoc(ClassName()+"::Sample_Binned (polymorphic)");
//        }
//
//        // Wrapper for the above function to facilitate runtime polymorphism.
//        virtual void Sample_Binned(
//            Real * restrict bins_out,
//            const Int bin_count_,
//            Real * restrict moments_out,
//            const Int moment_count_,
//            const Real * restrict ranges,
//            const std::vector< std::unique_ptr<RandomVariableBase_T> > & F_list_,
//            const Int sample_count,
//            const Int thread_count = 1
//        ) const override
//        {
//            const size_t size = F_list_.size();
//
//            std::vector< std::unique_ptr<RandomVariable_T> > F_list__ (size);
//
//            for( size_t i = 0; i < size; ++i )
//            {
//                F_list__[i] = std::unique_ptr<RandomVariable_T>(
//                    dynamic_cast<RandomVariable_T *>(
//                        F_list_[i]->Clone().get()
//                    )
//                );
//
//                print(F_list__[i]->Tag());
//
//                if( F_list__[i] == nullptr )
//                {
//                    eprint(ClassName()+"::Sample_Binned: Failed to downcast random variable "+F_list_[i]->Tag()+". Aborting.");
//
//                    return;
//                }
//            }
//
//            Sample_Binned( bins_out, bin_count_, moments_out, moment_count_, ranges, F_list__, sample_count, thread_count
//            );
//        }
//
//        virtual void NormalizeBinnedSamples(
//            Real * restrict bins,
//            const Int bin_count,
//            Real * restrict moments,
//            const Int moment_count,
//            const Int fun_count
//        ) const override
//        {
//            ptic(ClassName()+"::NormalizeBinnedSamples");
//            for( Int i = 0; i < 3; ++i )
//            {
//                for( Int j = 0; j < fun_count; ++j )
//                {
//                    // Normalize bins and moments.
//
//                    Real * restrict const bins_i_j = &bins[ (i*fun_count+j)*bin_count ];
//
//                    Real * restrict const moments_i_j = &moments[ (i*fun_count+j)*moment_count ];
//
//                    // The field for zeroth moment is assumed to contain the total mass.
//                    Real factor = Real(1)/moments_i_j[0];
//
//
//                    for( Int k = 0; k < bin_count; ++k )
//                    {
//                        bins_i_j[k] *= factor;
//                    }
//
//                    for( Int k = 0; k < moment_count; ++k )
//                    {
//                        moments_i_j[k] *= factor;
//                    }
//                }
//            }
//            ptoc(ClassName()+"::NormalizeBinnedSamples");
//        }
//
//#if defined(PLCTOPOLOGY_H)
//
//        std::map<std::string, std::tuple<Real,Real,Real>> SampleHOMFLY(
//            const Int sample_count,
//            const Int thread_count = 1
//        ) const
//        {
//            ptic(ClassName()+"SampleHOMFLY");
//
//            std::map<std::string, std::tuple<Real,Real,Real>> map_global;
//
//            if( AmbDim != 3 )
//            {
//                eprint("SampleHOMFLY is only available in 3D.");
//                return map_global;
//            }
//
//            JobPointers<Int> job_ptr (sample_count, thread_count);
//
//            valprint( "dimension   ", AmbDim       );
//            valprint( "edge_count  ", edge_count   );
//            valprint( "sample_count", sample_count );
//            valprint( "thread_count", thread_count );
//
//            gsl_rng_env_setup();
//
//            Tensor1<unsigned long,Int> seeds ( thread_count );
//
//            std::random_device r;
//
//            for( Int thread = 0 ; thread < thread_count; ++ thread )
//            {
//                seeds[thread] = (std::numeric_limits<unsigned int>::max()+1) * r() + r();
//            }
//
//            #pragma omp parallel for num_threads( thread_count )
//            for( Int thread = 0; thread < thread_count; ++thread )
//            {
//                const Int repetitions = (job_ptr[thread+1] - job_ptr[thread]);
//
//                Sampler W( edge_count, settings );
//
//                W.ReadEdgeLengths( r );
//                W.ReadRho( rho );
//
//                std::map<std::string, std::tuple<Real,Real,Real>> map_loc;
//
//                bool openQ = false;
//
//                int color_count = 0;
//
//                int n = static_cast<int>(edge_count);
//
//                gsl_rng *rng = gsl_rng_alloc( gsl_rng_mt19937 );
//
//                gsl_rng_set( rng, seeds[thread] );
//
//                plCurve * Gamma = plc_new( 1, &n, &openQ, &color_count );
////
////                auto * p = &Gamma->cp[0].vt[0];
//
//                auto & p = W.SpaceCoordinates();
//
//                for( Int l = 0; l < repetitions; ++l )
//                {
////                    valprint("l",l);
//
//                    W.RandomizeInitialEdgeCoordinates();
//
//                    W.ComputeShiftVector();
//
//                    W.Optimize();
//
//                    W.ComputeSpaceCoordinates();
//
//                    W.ComputeEdgeSpaceSamplingWeight();
//
//                    W.ComputeEdgeQuotientSpaceSamplingCorrection();
//
//                    const Real K = W.EdgeSpaceSamplingWeight();
//
//                    const Real K_quot = W.EdgeQuotientSpaceSamplingWeight();
//
//                    auto * comp = &Gamma->cp[0];
//
//                    for( int k = 0; k < edge_count; ++k )
//                    {
//                        auto * v = &comp->vt[k].c[0];
//
//                        for( int i = 0; i < 3; ++i )
//                        {
//                            v[i] = p(k,i);
//                        }
//                    }
//
//                    plc_fix_wrap( Gamma );
//
//
////                    char * polynomial = (char *) malloc( sizeof(char) * 3);
////
////                    polynomial[0] = 'a';
////                    polynomial[1] = 'b';
////                    polynomial[2] = 'c';
//
//                    char * polynomial = plc_homfly( rng, Gamma );
//
//                    std::string s("");
//
//                    if( polynomial != nullptr )
//                    {
//                        s.append(polynomial);
//                        free( polynomial );
//                    }
//                    else
//                    {
//                        s.append("FAILED");
//                    }
//
////                    s << polynomial;
////
////                    std::string str = s.str();
////                    print(s);
//
//                    map_loc.try_emplace(s, std::tie(zero,zero,zero) );
//
//                    auto & tuple = map_loc[s];
//
//                    std::get<0>(tuple) += one;
//                    std::get<1>(tuple) += K;
//                    std::get<2>(tuple) += K_quot;
//
//
//                }
//
//                #pragma omp critical
//                {
//                    for ( auto const & [key, val] : map_loc )
//                    {
//                        map_global.try_emplace( key, std::tie(zero,zero,zero) );
//
//                        auto & from = val;
//                        auto & to   = map_global[key];
//
//                        std::get<0>(to) += std::get<0>(from);
//                        std::get<1>(to) += std::get<1>(from);
//                        std::get<2>(to) += std::get<2>(from);
//                    }
//                }
//
//                gsl_rng_free( rng );
//
//                plc_free( Gamma );
//            }
//
//
//            ptoc(ClassName()+"::SampleHOMFLY");
//
//            return map_global;
//        }
//
//#endif
//
//    public:
//
//        void RandomizeInitialEdgeCoordinates() override
//        {
//            for( Int k = 0; k < edge_count; ++k )
//            {
//                Real r2 = static_cast<Real>(0);
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    const Real z = normal_dist( random_engine );
//
//                    x(k,i) = z;
//
//                    r2 += z * z;
//                }
//
//                Real r_inv = one/std::sqrt(r2);
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    x(k,i) *= r_inv;
//                }
//            }
//        }
//
//
//        void RandomSphericalPoints(
//            Real * restrict x_out,
//            const Int sample_count,
//            const Int thread_count_ = 1
//        ) const override
//        {
//            const Int thread_count = (thread_count_<=0) ? omp_get_num_threads() : thread_count_;
//
//            if( thread_count == 1 )
//            {
//                for( Int k = 0; k < edge_count; ++k )
//                {
//                    Real r2 = static_cast<Real>(0);
//
//                    Real v [AmbDim];
//
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        v[i] = normal_dist( random_engine );
//
//                        r2 += v[i] * v[i];
//                    }
//
//                    Real r_inv = one/std::sqrt(r2);
//
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        x_out[ AmbDim * k + i ] = v[i] * r_inv;
//                    }
//                }
//            }
//            else
//            {
//                JobPointers<Int> job_ptr ( sample_count, thread_count );
//
//                #pragma omp parallel for num_threads( thread_count )
//                for( Int thread = 0; thread < thread_count; ++thread )
//                {
//                    Real v [AmbDim];
//
//                    std::random_device r;
//
//                    std::seed_seq seed { r(), r(), r(), r() };
//
//                    std::mt19937_64 random_engine_loc { seed };
//
//                    std::normal_distribution<Real> dist { static_cast<Real>(0),static_cast<Real>(1) };
//
//                    const Int l_begin = job_ptr[thread];
//                    const Int l_end   = job_ptr[thread+1];
//
//                    for( Int l = l_begin; l < l_end; ++l )
//                    {
//                        Real * restrict x = &x_out[AmbDim * edge_count * l];
//
//                        for( Int k = 0; k < edge_count; ++k )
//                        {
//                            Real r2 = static_cast<Real>(0);
//
//                            for( Int i = 0; i < AmbDim; ++i )
//                            {
//                                v[i] = dist( random_engine_loc );
//
//                                r2 += v[i] * v[i];
//                            }
//
//                            Real r_inv = one/std::sqrt(r2);
//
//                            for( Int i = 0; i < AmbDim; ++i )
//                            {
//                                x[ AmbDim * k + i ] = v[i] * r_inv;
//                            }
//                        }
//                    }
//                }
//            }
//        }
//
//    protected:
//
//        ShiftMap<AmbDim,Real,Int> & Shifter()
//        {
//            return S;
//        }
//
//    protected:
//
//        Real tanhc( const Real t ) const
//        {
//            // Computes tanh(t)/t in a stable way by using a Padé approximation around t = 0.
//            constexpr Real a0 = static_cast<Real>(1);
//            constexpr Real a1 = static_cast<Real>(7)/static_cast<Real>(51);
//            constexpr Real a2 = static_cast<Real>(1)/static_cast<Real>(255);
//            constexpr Real a3 = static_cast<Real>(2)/static_cast<Real>(69615);
//            constexpr Real a4 = static_cast<Real>(1)/static_cast<Real>(34459425);
//
//            constexpr Real b0 = static_cast<Real>(1);
//            constexpr Real b1 = static_cast<Real>(8)/static_cast<Real>(17);
//            constexpr Real b2 = static_cast<Real>(7)/static_cast<Real>(255);
//            constexpr Real b3 = static_cast<Real>(4)/static_cast<Real>(9945);
//            constexpr Real b4 = static_cast<Real>(1)/static_cast<Real>(765765);
//
//            const Real t2 = t * t;
//
//            const Real result = ( t2 <= one ) ? (
//                a0 + t2 * (a1 + t2 * (a2 + t2 * (a3 + t2 * a4)))
//            )/(
//                b0 + t2 * (b1 + t2 * (b2 + t2 * (b3 + t2 * b4)))
//            )
//            : ( t2 <= static_cast<Real>(7) ) ? std::tanh(t)/t : one/std::abs(t);
//
//            return result;
//        }
//
//    public:
//
//        Int AmbientDimension() const override
//        {
//            return AmbDim;
//        }
//
//        virtual std::string ClassName() const override
//        {
//            return "Sampler<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+","+">";
//        }
//    };
//
//} // namespace CycleSampler
//
//#undef x
//#undef y
//#undef p
