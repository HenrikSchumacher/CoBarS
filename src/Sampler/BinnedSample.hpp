public:

    virtual void BinnedSample(
        mptr<Real> bins_inout,
        const Int bin_count_,
        mptr<Real> moments_inout,
        const Int moment_count_,
        cptr<Real> ranges,
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        // This function does the sampling, but computes moments and binning on the fly, so that the sampled data can be discarded immediately.
        
        // moments: A 3D-array of size 3 x fun_count x bin_count. Entry moments(i,j,k) will store the sampled weighted k-th moment of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i (see above).
        
        // ranges: Specify the range for binning: For j-th function in F_list, the range from ranges(j,0) to ranges(j,1) will be divided into bin_count bins. The user is supposed to provide meaningful ranges. Some rough guess might be obtained by calling the random variables on the prepared Sampler S.
        
        ptic(ClassName()+"::BinnedSample");
        
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
        
        
//        Tensor3<Real,Int> bins_global   ( bins_out,    3, fun_count, bin_count    );
//        Tensor3<Real,Int> moments_global( moments_out, 3, fun_count, moment_count );
        Tensor1<Real,Int> factor        (                 fun_count               );
        
        print("Sampling (binned) the following random variables:");
        for( Int i = 0; i < fun_count; ++ i )
        {
            const Size_T i_ = static_cast<Size_T>(i);
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
                
                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                
                S.LoadRandomVariables( F_list_ );
                
                Tensor3<Real,Int> bins_local   ( 3, fun_count, bin_count,    zero );
                Tensor3<Real,Int> moments_local( 3, fun_count, moment_count, zero );
                
                for( Int k = 0; k < repetitions; ++k )
                {
                    S.RandomizeInitialEdgeCoordinates();

                    S.ComputeConformalClosure();

                    S.ComputeSpaceCoordinates();
                    
                    S.ComputeEdgeSpaceSamplingWeight();
                    
                    S.ComputeEdgeQuotientSpaceSamplingCorrection();
                    
                    const Real K = S.EdgeSpaceSamplingWeight();

                    const Real K_quot = S.EdgeQuotientSpaceSamplingWeight();

                    for( Int i = 0; i < fun_count; ++i )
                    {
                        const Real val = S.EvaluateRandomVariable(i);

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
                        bins_local.data(), bins_inout, 3 * fun_count * bin_count
                    );
                    
                    add_to_buffer<VarSize,Sequential>(
                        moments_local.data(), moments_inout, 3 * fun_count * moment_count
                    );
                }
                
                Time stop = Clock::now();
             
                logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Duration(start, stop) ) + "." );
                
            },
            thread_count
        );
        
        ptoc(ClassName()+"::BinnedSample");
    }

