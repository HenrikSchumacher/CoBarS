public:
    
    virtual void Sample(
        Real * restrict const sampled_values,
        Real * restrict const K_edge_space,
        Real * restrict const K_quot_space,
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list,
        const Int sample_count,
        const Int  thread_count = 1
    ) const override
    {   
        TOOLS_PTIC(ClassName()+"::Sample (batch)");

        if( K_edge_space != nullptr )
        {
            sample_1<true>(
                sampled_values, K_edge_space, K_quot_space,
                F_list, sample_count, thread_count
            );
        }
        else
        {
            sample_1<false>(
                sampled_values, K_edge_space, K_quot_space,
                F_list, sample_count, thread_count
            );
        }
        
        TOOLS_PTOC(ClassName()+"::Sample (batch)");
    }


    virtual void Sample(
        Real * restrict const sampled_values,
        Real * restrict const K_edge_space,
        Real * restrict const K_quot_space,
        std::shared_ptr<RandomVariable_T> & F,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        // This function creates samples for the random variable F and records the sampling weights, so that this weighted data can be processed elsewhere.
        
        // The generated polygons are discarded immediately after evaluating the random variables on them.
        
        // sampled_values is expected to be an array of size at least sample_count;
        // K_edge_space is expected to be an array of size at least sample_count -- or a nullptr.
        // K_quot_space is expected to be an array of size at least sample_count  -- or a nullptr.
        
        TOOLS_PTIC(ClassName()+"::Sample");
        
        std::vector< std::shared_ptr<RandomVariable_T>> F_list;
        
        F_list.push_back( F->Clone() );

        Sample(
            sampled_values, K_edge_space, K_quot_space,
            F_list, sample_count, thread_count
        );
        
        TOOLS_PTOC(ClassName()+"::Sample");
    }


private:

    template<bool edge_space_flag>
    void sample_1(
        mptr<Real> sampled_values,
        mptr<Real> K_edge_space,
        mptr<Real> K_quot_space,
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list,
        const Int sample_count,
        const Int thread_count = 1
    ) const
    {
        if( K_quot_space != nullptr )
        {
            sample_2<edge_space_flag,true>(
                sampled_values, K_edge_space, K_quot_space,
                F_list, sample_count, thread_count
            );
        }
        else
        {
            sample_2<edge_space_flag,false>(
                sampled_values, K_edge_space, K_quot_space,
                F_list, sample_count, thread_count
            );
        }
    }


    template<bool edge_space_flag, bool quotient_space_flag>
    void sample_2(
        mptr<Real> sampled_values,
        mptr<Real> K_edge_space,
        mptr<Real> K_quot_space,
        const std::vector< std::shared_ptr<RandomVariable_T> > & F_list,
        const Int sample_count,
        const Int thread_count = 1
    ) const
    {
        const Int fun_count = static_cast<Int>(F_list.size());
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );

                // For every thread create a copy of the current Sampler object.
                Sampler S ( EdgeLengths().data(), Rho().data(), EdgeCount(), Settings() );
                
                S.LoadRandomVariables( F_list );
                
                for( Int k = k_begin; k < k_end; ++k )
                {
                    S.RandomizeInitialEdgeVectors();
                    
                    S.ComputeConformalClosure();
                    
                    if constexpr ( edge_space_flag )
                    {
                        K_edge_space[k] = S.EdgeSpaceSamplingWeight();
                    }
                    
                    if constexpr ( quotient_space_flag )
                    {
                        K_quot_space[k] = S.EdgeQuotientSpaceSamplingWeight();
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

