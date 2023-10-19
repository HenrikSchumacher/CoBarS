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
            sample_1<true>(
                sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                F_list_, sample_count, thread_count
            );
        }
        else
        {
            sample_1<false>(
                sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                F_list_, sample_count, thread_count
            );
        }
        
        ptoc(ClassName()+"::Sample (batch)");
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
        // This function creates samples for the random variable F and records the sampling weights, so that this weighted data can be processed elsewhere.
        
        // The generated polygons are discarded immediately after evaluating the random variables on them.
        
        // sampled_values is expected to be an array of size at least sample_count;
        // edge_space_sampling_weights is expected to be an array of size at least sample_count -- or a nullptr.
        // edge_quotient_space_sampling_weights is expected to be an array of size at least sample_count  -- or a nullptr.
        
        ptic(ClassName()+"::Sample");
        
        std::vector< std::shared_ptr<RandomVariable_T>> F_list_;
        
        F_list_.push_back( F_->Clone() );

        Sample(
            sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
            F_list_, sample_count, thread_count
        );
        
        ptoc(ClassName()+"::Sample");
    }


private:

    template<bool edge_space_flag>
    void sample_1(
        mut<Real> sampled_values,
        mut<Real> edge_space_sampling_weights,
        mut<Real> edge_quotient_space_sampling_weights,
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
        const Int sample_count,
        const Int thread_count = 1
    ) const
    {
        if( edge_quotient_space_sampling_weights != nullptr )
        {
            sample_2<edge_space_flag,true>(
                sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                F_list_, sample_count, thread_count
            );
        }
        else
        {
            sample_2<edge_space_flag,false>(
                sampled_values, edge_space_sampling_weights, edge_quotient_space_sampling_weights,
                F_list_, sample_count, thread_count
            );
        }
    }


    template<bool edge_space_flag, bool quotient_space_flag>
    void sample_2(
        mut<Real> sampled_values,
        mut<Real> edge_space_sampling_weights,
        mut<Real> edge_quotient_space_sampling_weights,
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
        const Int sample_count,
        const Int thread_count = 1
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
                
                S.LoadRandomVariables( F_list_ );
                
                for( Int k = k_begin; k < k_end; ++k )
                {
                    S.RandomizeInitialEdgeCoordinates();
                    
                    S.ComputeConformalClosure();
                    
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
                        edge_quotient_space_sampling_weights[k] = S.EdgeQuotientSpaceSamplingWeight();
                    }
                    
                    for( Int i = 0; i < fun_count; ++i )
                    {
                        sampled_values[k * fun_count + i] = S.EvaluateRandomVariable(i);
                    }
                }
            },
            thread_count
        );
    }

