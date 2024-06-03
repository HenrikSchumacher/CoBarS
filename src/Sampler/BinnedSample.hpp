public:

    virtual void BinnedSample(
              Real * restrict const bins,   const Int bin_count,
              Real * restrict const moms,   const Int mom_count,
        const Real * restrict const ranges,
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        // This function does the sampling, but computes moments and binning on the fly, so that the sampled data can be discarded immediately.
        
        // moments: A 3D-array of size 3 x fun_count x bin_count. Entry moments(i,j,k) will store the sampled weighted k-th moment of the j-th random variable from the list F_list -- with respect to the weights corresponding to the value of i (see above).
        
        // ranges: Specify the range for binning: For j-th function in F_list, the range from ranges(j,0) to ranges(j,1) will be divided into bin_count bins. The user is supposed to provide meaningful ranges. Some rough guess might be obtained by calling the random variables on the prepared Sampler S.
        
        ptic(ClassName()+"::BinnedSample");
        
        const Int f_count = static_cast<Int>(F_list.size());
        
        const Int m_count = std::max( static_cast<Int>(3), mom_count );
        
        const Int b_count = std::max( bin_count, static_cast<Int>(1) );
        
        valprint( "dimension   ", AmbDim       );
        valprint( "edge_count  ", edge_count_  );
        valprint( "sample_count", sample_count );
        valprint( "fun_count   ", f_count      );
        valprint( "bin_count   ", b_count      );
        valprint( "moment_count", m_count      );
        valprint( "thread_count", thread_count );
        
        
        Tensor1<Real,Int> factor ( f_count );
        
        print("Sampling (binned) the following random variables:");
        for( Int i = 0; i < f_count; ++ i )
        {
            const Size_T i_ = static_cast<Size_T>(i);
            factor(i) = static_cast<Real>(bin_count) / ( ranges[2*i+1] - ranges[2*i+0] );

            print("    " + F_list[i_]->Tag());
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
                
                Sampler S ( EdgeLengths().data(), Rho().data(), EdgeCount(), Settings() );
                
                S.LoadRandomVariables( F_list );

                Tensor3<Real,Int> bins_local( 3, f_count, b_count, zero );
                Tensor3<Real,Int> moms_local( 3, f_count, m_count, zero );

                for( Int k = 0; k < repetitions; ++k )
                {
                    S.RandomizeInitialEdgeVectors();

                    S.ComputeConformalClosure();
                    
                    const Real K = S.EdgeSpaceSamplingWeight();

                    const Real K_quot = S.EdgeQuotientSpaceSamplingWeight();

                    for( Int i = 0; i < f_count; ++i )
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

                        moms_local(0,i,0) += values[0];
                        moms_local(1,i,0) += values[1];
                        moms_local(2,i,0) += values[2];

                        for( Int j = 1; j < m_count; ++j )
                        {
                            values[0] *= val;
                            values[1] *= val;
                            values[2] *= val;
                            moms_local(0,i,j) += values[0];
                            moms_local(1,i,j) += values[1];
                            moms_local(2,i,j) += values[2];
                        }
                    }
                }
                
                {
                    const std::lock_guard<std::mutex> lock ( mutex );
                    
                    add_to_buffer<VarSize,Sequential>(
                        bins_local.data(), bins, 3 * f_count * b_count
                    );
                    
                    add_to_buffer<VarSize,Sequential>(
                        moms_local.data(), moms, 3 * f_count * m_count
                    );
                }
                
                Time stop = Clock::now();
             
                logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Duration(start, stop) ) + "." );
                
            },
            thread_count
        );
        
        ptoc(ClassName()+"::BinnedSample");
    }

