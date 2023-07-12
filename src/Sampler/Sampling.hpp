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

