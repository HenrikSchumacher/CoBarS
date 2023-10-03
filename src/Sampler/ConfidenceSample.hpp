public:

    // TODO: Add user-defined max_sample_count

    // TODO: Use long double for moment containers?


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
        
        moments = Tensor2<Real,Int> ( 3, fun_count + 1, zero );
        
        // moments[0] - 1-st moments
        // moments[1] - 2-nd moments
        // moments[2] - moments for covariance
        
        // Prepare samplers.
        
        std::vector<Sampler> samplers (thread_count);
        
        for( Int thread = 0; thread < thread_count; ++thread )
        {


        }
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                                
                S.LoadRandomVariables( F_list_ );
                
                S.moments = Tensor2<Real,Int>( 3, fun_count + 1, zero );
                
                samplers[thread] = std::move(S);
                
            },
            thread_count
        );
        
        
        ptoc("Preparation");
        
        std::mutex mutex;
        
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

                        S.ComputeShiftVector();

                        S.Optimize();

                        S.ComputeSpaceCoordinates();
                        
                        S.ComputeEdgeSpaceSamplingWeight();
                        
                        Real weight = 0;
                        
                        if constexpr ( quotient_space_Q )
                        {
                            S.ComputeEdgeQuotientSpaceSamplingCorrection();
                            
                            weight = S.EdgeQuotientSpaceSamplingWeight();
                        }
                        else
                        {
                            weight = S.EdgeSpaceSamplingWeight();
                        }
                        
                        Real weight_squared = weight * weight;
                        
                        for( Int i = 0; i < fun_count; ++i )
                        {
                            const Real value = S.EvaluateRandomVariable(i);

                            const Real value_squared = value * value;
                            
                            S.moments[0][i] += weight * value;
                            S.moments[1][i] += weight * value_squared;
                            S.moments[2][i] += weight_squared * value;
                        }
                        
                        S.moments[0][fun_count] += weight;
                        S.moments[1][fun_count] += weight_squared;
                    }
                    
                    {
                        const std::lock_guard<std::mutex> lock ( mutex );
                        
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
            
            const Real mean_K = Frac( moments[0][fun_count], N );
            
            const Real var_K  = Bessel_corr * ( Frac( moments[1][fun_count], N ) - mean_K * mean_K );
            
            const Real mean_Y = mean_K;
            
            const Real var_Y  = Frac( var_K, N );
            
//            dump( mean_Y );
//            dump( std::sqrt(var_Y)  );
            
            // Check Geary condition
            if( mean_Y < Scalar::Three<Real> * std::sqrt( var_Y ) )
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
                    const Real mean_F = Frac( moments[0][i], N );
                    
                    const Real var_F  = Bessel_corr * ( Frac( moments[1][i], N ) - mean_F * mean_F );
                    
                    const Real cov_FK = Bessel_corr * ( Frac( moments[2][i], N ) - mean_F * mean_K );
                    
                    const Real T = moments[0][i] / moments[0][fun_count];

                    
                    const Real mean_X = mean_F;
                    
                    const Real var_X  = Frac( var_F, N );
                    
                    const Real cov_XY = Frac( cov_FK, N );
                    
//                    dump( mean_X );
//                    dump( std::sqrt(var_X)  );
//                    dump( cov_XY  / ( std::sqrt(var_X * var_Y) ));
                    
                    const GearyTransform G ( mean_X, mean_Y, var_X, cov_XY, var_Y );
                                    
                    // t - G( T ) is a standard Gaussian.
                    // See Geary - The Frequency Distribution of the Quotient of Two Normal Variates (1930)
                    // https://www.jstor.org/stable/2342070
                    
                    const Real z_lo = G( T - radii[i] );
                    const Real z_hi = G( T + radii[i] );
                    
                    // Check for sufficient confidence.
                    
//                    dump( T );
//                    dump( z_lo );
//                    dump( z_hi );
                    
                    const Real current_confidence = N_CDF(z_hi) - N_CDF(z_lo);
                    
                    valprint( "current confidence of " + F_list_[i]->Tag(), current_confidence );
                                        
                    completed = completed && ( current_confidence > confidence );
                }
            }
            
        }
        
        ptoc("Sampling");
        
        if( !completed )
        {
            // TODO: print warnings.
        }
        
        valprint( "N", N );
        
        ptic("Postprocessing");
        
        const Real Bessel_corr = Frac<Real>( N, N-1);
        
        const Real mean_K = Frac( moments[0][fun_count], N );
        
        const Real var_K  = Bessel_corr * ( Frac( moments[1][fun_count], N ) - mean_K * mean_K );
        
        const Real mean_Y = mean_K;
        
        const Real var_Y  = Frac( var_K, N );
        
        for( Int i = 0; i < fun_count; ++i )
        {
            print( F_list_[i]->Tag() );
            
            const Real mean_F = Frac( moments[0][i], N );
            
            const Real var_F  = Bessel_corr * ( Frac( moments[1][i], N ) - mean_F * mean_F );
            
            const Real cov_FK = Bessel_corr * ( Frac( moments[2][i], N ) - mean_F * mean_K );
            
            const Real mean_X = mean_F;
            
            const Real var_X  = Frac( var_F, N );
            
            const Real cov_XY = Frac( cov_FK, N );
            
            const GearyTransform G ( mean_X, mean_Y, var_X, cov_XY, var_Y );
            
            const Real T = moments[0][i] / moments[0][fun_count];
            
            
            // Newton's method to find confidence interval
            // Initial guess.
            
            Real error_F = radii[i];
            
              // TODO: Add line search!
            {
                Real T_lo = T - error_F;
                Real T_hi = T + error_F;
                
                Real z_lo = G( T_lo );
                Real z_hi = G( T_hi );
                
                Real failure = N_CDF(z_hi) - N_CDF(z_lo) - confidence;
                
                dump(failure);
                
                Int Newton_iter = 0;
                Int max_Newton_iter = 10;
                
                while( ( std::abs(failure) >= 0.0001 ) && ( Newton_iter < max_Newton_iter ) )
                {
                    ++Newton_iter;
                    const Real Dfailure = N_PDF(z_hi) * G.D(T_hi) + N_PDF(z_lo) * G.D(T_lo);
                    
                    dump(Dfailure);
                    
                    error_F -= failure / Dfailure;
                    
                    T_lo = T - error_F;
                    T_hi = T + error_F;
                    
                    z_lo = G( T_lo );
                    z_hi = G( T_hi );
                    
                    failure = N_CDF(z_hi) - N_CDF(z_lo) - confidence;
                    
                    dump(failure);
                }
                
                // TODO: Better warning message.
                if( Newton_iter >= max_Newton_iter )
                {
                    dump(Newton_iter);
                    
                    wprint( "max_Newton_iter" );
                }
            }
            
            means[i]  = mean_F / mean_K;
            errors[i] = error_F;
            
//            dump( means[i] );
//            dump( errors[i] );
        }
        
        ptoc("Postprocessing");
        
        moments = Tensor2<Real,Int>();
        
        ptoc(ClassName()+"::ConfidenceSample");
        
        return N;
    }




