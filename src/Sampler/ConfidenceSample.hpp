public:
   
    virtual Int ConfidenceSample(
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
        mptr<Real> means,
        mptr<Real> errors,
        cptr<Real> radii,  // desired radii of the confidence intervals
        const Int  max_sample_count,
        const bool quotient_space_Q,
        const Int thread_count = 1,
        const Real confidence = 0.95,
        const Int chunk_size = 1000000
    ) const override
    {
        // means, errors and radii are expected to be allocated arrays sufficiently large to hold at least F_list_.size() numbers.
        
        if ( quotient_space_Q )
        {
            return confidenceSample<true>(
                F_list_,
                means, errors, radii,
                max_sample_count, thread_count, confidence, chunk_size
            );
        }
        else
        {
            return confidenceSample<false>(
                F_list_,
                means, errors, radii,
                max_sample_count, thread_count, confidence, chunk_size
            );
        }
    }


private:

    template<bool quotient_space_Q>
    Int confidenceSample(
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
        mptr<Real> means,
        mptr<Real> errors,
        cptr<Real> radii,  // desired radii of the confidence intervals
        const Int  max_sample_count,
        const Int  thread_count = 1,
        const Real confidence = 0.95,
        const Int  chunk_size = 1000000
    ) const
    {
        // Samples the random variables in F_list_ until the radii of the confidence intervals of each function are lower or equal to the prescribed radii (or until max_sample_count samples have been drawn, whatever happens first).
        
        ptic(ClassName()+"::ConfidenceSample");
        
        if( (confidence < zero) )
        {
            eprint(ClassName()+"::ConfidenceSample: confidence level " + ToString(confidence) + " is smaller than zero. Aborting." );
            
            return 0;
        }
        
        if( (confidence > one) )
        {
            eprint(ClassName()+"::ConfidenceSample: confidence level " + ToString(confidence) + " is greater than 1. Aborting." );
            
            return 0;
        }
        
        if( (confidence > small_one) )
        {
            wprint(ClassName()+"::ConfidenceSample: confidence level " + ToString(confidence) + " is too close to 1. Computing max_sample_count samples" );
//
//            return 0;
        }
        
        ptic("Preparation");
        
        const Int fun_count = static_cast<Int>(F_list_.size());
                
                                                            
        valprint( "dimension       ", AmbDim            );
        valprint( "edge_count      ", edge_count        );
        valprint( "fun_count       ", fun_count         );
        valprint( "thread_count    ", thread_count      );
        valprint( "confidence level", confidence        );
        
        print("ConfidenceSample works on the following random variables:");
        
        for( RandomVariable_Ptr f : F_list_ )
        {
            print("    " + f->Tag());
        }
        
        Int N = 0;
        
        moments.Resize( 3, fun_count + 1 );
        moments.SetZero();
        
        // moments[0] - 1-st moments
        // moments[1] - 2-nd moments
        // moments[2] - moments for covariance
        
        // Prepare samplers.
        
        std::vector<Sampler> samplers (thread_count);
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                                
                S.LoadRandomVariables( F_list_ );
                
                S.moments.Resize( 3, fun_count + 1 );
                
                samplers[thread] = std::move(S);
                
            },
            thread_count
        );
        
        
        ptoc("Preparation");
        
        std::mutex moment_mutex;
        
        bool completed = false;
        
        ptic("Sampling");
        
        while( !completed )
        {
            ParallelDo(
                [&,this]( const Int thread )
                {
                    Time start = Clock::now();
                    
                    const Int k_begin = JobPointer( chunk_size, thread_count, thread     );
                    const Int k_end   = JobPointer( chunk_size, thread_count, thread + 1 );
                    
                    const Int repetitions = k_end - k_begin;
                    
                    Sampler & S = samplers[thread];
                    
                    S.moments.SetZero();
                    
                    for( Int k = 0; k < repetitions; ++k )
                    {
                        S.RandomizeInitialEdgeCoordinates();

                        S.ComputeConformalClosure();

                        S.ComputeSpaceCoordinates();
                        
                        S.ComputeEdgeSpaceSamplingWeight();
                        
                        Real K = 0;
                        
                        if constexpr ( quotient_space_Q )
                        {
                            S.ComputeEdgeQuotientSpaceSamplingCorrection();
                            
                            K = S.EdgeQuotientSpaceSamplingWeight();
                        }
                        else
                        {
                            K = S.EdgeSpaceSamplingWeight();
                        }
                        
                        for( Int i = 0; i < fun_count; ++i )
                        {
                            const Real KF = K * S.EvaluateRandomVariable(i);
                            
                            S.moments[0][i] += KF;
                            S.moments[1][i] += KF * KF;
                            S.moments[2][i] += KF * K ;

                        }
                        
                        S.moments[0][fun_count] += K;
                        S.moments[1][fun_count] += K * K;
                    }
                    
                    {
                        const std::lock_guard<std::mutex> lock ( moment_mutex );
                        
                        add_to_buffer(
                            S.moments.data(), moments.data(), 3 * (fun_count+1)
                        );
                    }
                    
                    Time stop = Clock::now();
                 
                    logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Duration(start, stop) ) + "." );
                    
                },
                thread_count
            );
            
            N += chunk_size;
            
            dump( N );
            
            if( N > max_sample_count )
            {
                wprint(ClassName()+"::ConfidenceSample: Maximal number of samples reached. Sampling aborted after " + ToString(N) + " samples.");
                break;
            }
            
            const Real Bessel_corr = Frac<Real>( N, N-1 );
            
            const Real mean_K = Frac<Real>( moments[0][fun_count], N );
            
            const Real var_K  = Bessel_corr * ( Frac<Real>( moments[1][fun_count], N ) - mean_K * mean_K );
            
            const Real mean_Y = mean_K;
            
            const Real var_Y  = Frac<Real>( var_K, N );
            
//            dump( mean_Y );
//            dump( std::sqrt(var_Y)  );
            
            const Real Geary_factor = mean_Y / std::sqrt( var_Y );
            
            dump( Geary_factor );
            
            // Check Geary condition
            if( Geary_factor < Scalar::Three<Real> )
            {
                wprint("Geary condition failed.");
                
                completed = false;
            }
            else
            {
                
                completed = true;
                
                // Check stopping criterion for each random variable.
                                
                for( Int i = 0; i < fun_count; ++i )
                {
                    const Real mean_KF  = Frac<Real>( moments[0][i], N );
                    
                    const Real var_KF   = Bessel_corr * ( Frac<Real>( moments[1][i], N ) - mean_KF * mean_KF );
                    
                    const Real cov_KF_K = Bessel_corr * ( Frac<Real>( moments[2][i], N ) - mean_KF * mean_K  );
                    
                    const Real mean_X = mean_KF;
                    
                    const Real var_X  = Frac<Real>( var_KF, N );
                    
                    const Real cov_XY = Frac<Real>( cov_KF_K, N );
                    
                    const Real T = mean_X / mean_Y;
                    
                    // TODO: Boundary of checked world.
                    
                    const GearyTransform<Real> G ( mean_X, mean_Y, var_X, cov_XY, var_Y );
                                    
                    // t - G( T ) is a standard Gaussian.
                    // See Geary - The Frequency Distribution of the Quotient of Two Normal Variates (1930)
                    // https://www.jstor.org/stable/2342070
                    
                    const Real z_lo = G( T - radii[i] );
                    const Real z_hi = G( T + radii[i] );
                    
                    // Check for sufficient confidence.
                    
                    const Real current_confidence = N_CDF(z_hi) - N_CDF(z_lo);

                    
//                    valprint( "current confidence of  " + F_list_[i]->Tag(), current_confidence );
                    
                    print( "  Current estimate of " + F_list_[i]->Tag() + " = " +  ToString(T) + " +/- " + ToString(radii[i]) + " with confidence = " + ToString(current_confidence) + "." );
                                        
                    completed = completed && ( current_confidence > confidence );
                }
            }
            
        }
        
        ptoc("Sampling");
        
        
        valprint( "N", N );
        
        ptic("Postprocessing");
        
        const Real Bessel_corr = Frac<Real>( N, N-1 );
        
        const Real mean_K = Frac( moments[0][fun_count], N );
        
        const Real var_K  = Bessel_corr * ( Frac( moments[1][fun_count], N ) - mean_K * mean_K );
        
        const Real mean_Y = mean_K;
        
        const Real var_Y  = Frac( var_K, N );
        
        
        for( Int i = 0; i < fun_count; ++i )
        {
            print( F_list_[i]->Tag() );
            
            const Real mean_KF  = Frac<Real>( moments[0][i], N );
            
            const Real var_KF   = Bessel_corr * ( Frac<Real>( moments[1][i], N ) - mean_KF * mean_KF );
            
            const Real cov_KF_K = Bessel_corr * ( Frac<Real>( moments[2][i], N ) - mean_KF * mean_K  );
            
            const Real mean_X = mean_KF;

            const Real var_X  = Frac<Real>( var_KF, N );
            
            const Real cov_XY = Frac<Real>( cov_KF_K, N );
            
            const GearyTransform<Real> G ( mean_X, mean_Y, var_X, cov_XY, var_Y );
            
            const Real T = mean_X / mean_Y;
            
            // Interval bisection to find the actual confidence radius.
            
            Real a = 0;
            Real b = radii[i];
            
            Real P_a = 0;
            Real P_b = N_CDF( G( T + b ) ) - N_CDF( G( T - b ) );
            
            
            // Extend the search interval to make sure that the actual confidence radius lies within [a,b)
            while( P_b <= confidence )
            {
                a = b;
                
                b *= Scalar::Two<Real>;
                
                P_b = N_CDF( G( T + b ) ) - N_CDF( G( T - b ) );
            }
            
            while( b - a > 0.001 * b )
            {
                const Real c = Scalar::Half<Real> * (a + b);
                
                const Real P_c = N_CDF( G( T + c ) ) - N_CDF( G( T - c ) );
                
                if( P_c > confidence )
                {
                    b   = c;
                    P_b = P_c;
                }
                else
                {
                    a   = c;
                    P_a = P_c;
                }
            }
            
            means[i]  = mean_KF / mean_K;
            errors[i] = Scalar::Half<Real> * b;
        }
        
        ptoc("Postprocessing");
        
        moments.Resize(0,0);
        
        ptoc(ClassName()+"::ConfidenceSample");
        
        return N;
    }




