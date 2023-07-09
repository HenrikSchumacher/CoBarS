public:

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
        
        if constexpr ( zerofyfirstQ )
        {
            F.SetZero();
            DF.SetZero();
            
            for( Int k = 0; k < edge_count; ++k )
            {
                Vector_T y_k ( y, k );
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = r[k] * y_k[i];
                    
                    F[i] -= factor;
                    
                    for( Int j = i; j < AmbDim; ++j )
                    {
                        DF[i][j] -= factor * y_k[j];
                    }
                }
            }
        }
        else
        {
            // Filling F and DF with first summand...
            {
                Vector_T y_k ( y, 0 );
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = r[0] * y_k[i];
                    
                    F[i] = - factor;
                    
                    for( Int j = i; j < AmbDim; ++j )
                    {
                        DF[i][j] = - factor * y_k[j];
                    }
                }
            }
            
            // ... and adding-in the other summands.
            for( Int k = 1; k < edge_count; ++k )
            {
                Vector_T y_k ( y, k );
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = r[k] * y_k[i];
                    
                    F[i] -= factor;
                    
                    for( Int j = i; j < AmbDim; ++j )
                    {
                        DF[i][j] -= factor * y_k[j];
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
                
                Vector_T y_k;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    y_k[i] = (one_minus_ww * x_k[i] + wx2_minus_2 * w[i]) * denom;
                }
                
                y_k.Write( y, k );
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

                Sampler S( edge_count, Settings() );

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

