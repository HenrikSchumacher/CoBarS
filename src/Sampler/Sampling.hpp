public:

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

                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );

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

                // For every thread create a copy of the current Sampler object.
                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                
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

                // For every thread create a copy of the current Sampler object.
                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                
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
        // ranges: Specify the range for binning: For j-th function in F_list, the range from ranges(j,0) to ranges(j,1) will be devided into bin_count bins. The user is supposed to provide meaningful ranges. Some rough guess might be obtained by calling the random variables on the prepared Sampler C.
        
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
                
                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );
                
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
