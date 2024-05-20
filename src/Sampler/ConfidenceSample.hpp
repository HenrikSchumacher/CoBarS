public:
   
    virtual Int ConfidenceSample(
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
        mptr<Real> sample_means,
        mptr<Real> sample_variances,
        mptr<Real> errors,
        cptr<Real> radii,  // desired radii of the confidence intervals
        const Int  max_sample_count,
        const bool quotient_space_Q,
        const Int  thread_count = 1,
        const Real confidence = 0.95,
        const Int  chunk_size = 1000000,
        const bool relativeQ = false,
        const bool verboseQ = true
    ) const override
    {
        // means, errors and radii are expected to be allocated arrays sufficiently large to hold at least F_list_.size() numbers.
        
        if ( quotient_space_Q )
        {
            return confidenceSample<true>(
                F_list_,
                sample_means, sample_variances, errors, radii,
                max_sample_count, thread_count, confidence, chunk_size, relativeQ, verboseQ
            );
        }
        else
        {
            return confidenceSample<false>(
                F_list_,
                sample_means, sample_variances, errors, radii,
                max_sample_count, thread_count, confidence, chunk_size, relativeQ, verboseQ
            );
        }
    }


private:


    template<bool quotient_space_Q>
    Int confidenceSample(
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
        mptr<Real> sample_means,
        mptr<Real> sample_variances,
        mptr<Real> errors,
        cptr<Real> radii,  // desired radii of the confidence intervals
        const Int  max_sample_count,
        const Int  thread_count,
        const Real confidence,
        const Int  chunk_size,
        const bool relativeQ,
        const bool verboseQ
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
                
              
        if( verboseQ )
        {
            valprint( "dimension       ", AmbDim            );
            valprint( "edge_count      ", edge_count        );
            valprint( "fun_count       ", fun_count         );
            valprint( "thread_count    ", thread_count      );
            valprint( "confidence level", confidence        );
            
            if( relativeQ )
            {
                print("Using relative error measures.");
            }
            else
            {
                print("Using absolute error measures.");
            }
            
            print("ConfidenceSample works on the following random variables:");
            
            for( RandomVariable_Ptr f : F_list_ )
            {
                print("    " + f->Tag());
            }
        }
        
        Int N = 0;
        
        Real total_time = 0;
        
        moments.Resize( 4, fun_count + 1 );
        moments.SetZero();
        
        // moments[0] - 1-st moments
        // moments[1] - 2-nd moments
        // moments[2] - moments for covariance
        // moments[3] - moments for sample variance
        
        // Prepare samplers.
        
        std::vector<Sampler> samplers (thread_count);
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                                
                S.LoadRandomVariables( F_list_ );
                
                S.moments.Resize( 4, fun_count + 1 );
                
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
            Time start_time = Clock::now();
            
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
                        
                        if constexpr ( quotient_space_Q )
                        {
                            S.ComputeEdgeQuotientSpaceSamplingCorrection();
                        }
                        
                        Real K = 0;
                        
                        if constexpr ( quotient_space_Q )
                        {
                            K = S.EdgeQuotientSpaceSamplingWeight();
                        }
                        else
                        {
                            K = S.EdgeSpaceSamplingWeight();
                        }
                        
                        for( Int i = 0; i < fun_count; ++i )
                        {
                            const Real F_ = S.EvaluateRandomVariable(i);
                            const Real KF = K * F_;
                            
                            S.moments[0][i] += KF;
                            S.moments[1][i] += KF * KF;
                            S.moments[2][i] += KF * K ;
                            S.moments[3][i] += KF * F_;

                        }
                        
                        S.moments[0][fun_count] += K;
                        S.moments[1][fun_count] += K * K;
                    }
                    
                    {
                        const std::lock_guard<std::mutex> lock ( moment_mutex );
                        
                        add_to_buffer(
                            S.moments.data(), moments.data(), 4 * (fun_count+1)
                        );
                    }
                    
                    Time stop = Clock::now();
                 
                    logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Duration(start, stop) ) + "." );
                    
                },
                thread_count
            );
            
            Time stop_time = Clock::now();
            
            Real time = Tools::Duration(start_time,stop_time);

            total_time += time;
            
            N += chunk_size;
            
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
            
            const Real Geary_factor = mean_Y / std::sqrt( var_Y );
            
            // Check Geary condition
            if( Geary_factor < static_cast<Real>(3) )
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
                    
                    const Real absolute_radius = relativeQ ? radii[i] * T : radii[i];
                    
                    // TODO: Boundary of checked world.
                    
                    const GearyTransform<Real> G ( mean_X, mean_Y, var_X, cov_XY, var_Y );
                                    
                    // t - G( T ) is a standard Gaussian.
                    // See Geary - The Frequency Distribution of the Quotient of Two Normal Variates (1930)
                    // https://www.jstor.org/stable/2342070
                    
                    const Real z_lo = G( T - absolute_radius );
                    const Real z_hi = G( T + absolute_radius );
                    
                    // Check for sufficient confidence.
                    
                    const Real current_confidence = N_CDF(z_hi) - N_CDF(z_lo);

                    if( verboseQ )
                    {
                        valprint("  total_time ", total_time );
                        
                        print( "  Current estimate of " + F_list_[i]->Tag() + " = " +  ToString(T) + " +/- " + ToString(absolute_radius) + " with confidence = " + ToString(current_confidence) + "." );
                    }
                    
                    completed = completed && ( current_confidence > confidence );
                }
            }
            
        }
        
        ptoc("Sampling");
        
        ptic("Postprocessing");
        
        const Real Bessel_corr = Frac<Real>( N, N-1 );
        
        const Real mean_K = Frac( moments[0][fun_count], N );
        
        const Real var_K  = Bessel_corr * ( Frac( moments[1][fun_count], N ) - mean_K * mean_K );
        
        const Real mean_Y = mean_K;
        
        const Real var_Y  = Frac( var_K, N );
        
        
        for( Int i = 0; i < fun_count; ++i )
        {
            const Real mean_KF  = Frac<Real>( moments[0][i], N );
            
            const Real var_KF   = Bessel_corr * ( Frac<Real>( moments[1][i], N ) - mean_KF * mean_KF );
            
            const Real cov_KF_K = Bessel_corr * ( Frac<Real>( moments[2][i], N ) - mean_KF * mean_K  );
            
            const Real mean_X = mean_KF;

            const Real var_X  = Frac<Real>( var_KF, N );
            
            const Real cov_XY = Frac<Real>( cov_KF_K, N );
            
            const GearyTransform<Real> G ( mean_X, mean_Y, var_X, cov_XY, var_Y );
            
            const Real T = mean_X / mean_Y;

            const Real absolute_radius = relativeQ ? radii[i] * T : radii[i];
            
            // Interval bisection to find the actual confidence radius.
            
            auto P = [T,&G]( const Real b )
            {
                return N_CDF( G( T + b ) ) - N_CDF( G( T - b ) );
            };
            
            Real a = 0;
            Real b = absolute_radius;
            
            // Extend the search interval to make sure that the actual confidence radius lies within [a,b)
            while( P(b) <= confidence )
            {
                a = b;
                
                b *= static_cast<Real>(2);
            }
            
            const Real error = BisectionSearch<1>( std::move(P), a, b, confidence, 0.001 );
            
            sample_means[i]     = T;
            sample_variances[i] = Bessel_corr * ( Frac( moments[3][i], moments[0][fun_count] ) - T * T );
            errors[i]           = error;
        }
        
        
        ptoc("Postprocessing");
        
        moments.Resize(0,0);
        
        ptoc(ClassName()+"::ConfidenceSample");
        
        return N;
    }




