private:
    // Common implementation of
    //
    // - ComputeConformalClosures,
    // - CreateRandomClosedPolygons,
    //
    // - ComputeConformalCentralizations,
    // - CreateRandomCentralizedPointClouds,
    // - CreateRandomCentralizedPointClouds_Detailed.
    //
    // This function is mean to be called only from there.

    template<
        bool x_in_Q, bool x_out_Q,
        bool p_in_Q, bool p_out_Q, bool w_Q, bool y_Q, bool q_Q,
        bool edge_space_Q, bool quot_space_Q
    >
    void CreatePolygons(
        cptr<Real> x_in, mptr<Real> x_out,
        cptr<Real> p_in, mptr<Real> p_out,
        mptr<Real> w,
        mptr<Real> y,
        mptr<Real> q,
        mptr<Real> K_edge_space,
        mptr<Real> K_quot_space,
        const Int sample_count,
        const Int thread_count
    ) const
    {
        ParallelDo(
            [&,this]( const Int thread )
            {
                Time start = Clock::now();
                
                const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );

                Sampler S ( EdgeLengths().data(), Rho().data(), EdgeCount(), Settings() );

                for( Int k = k_begin; k < k_end; ++k )
                {
                    if constexpr ( p_in_Q || x_in_Q )
                    {
                        if constexpr ( p_in_Q )
                        {
                            S.ReadInitialVertexPositions(p_in,k);
                        }
                        else
                        {
                            S.ReadInitialEdgeVectors(x_in,k);
                        }
                    }
                    else
                    {
                        S.RandomizeInitialEdgeVectors();
                    }
                    
                    if constexpr ( w_Q || y_Q || q_Q || edge_space_Q || quot_space_Q )
                    {
                        S.computeConformalClosure<q_Q,quot_space_Q>();
                    }
                    
                    
                    if constexpr ( p_out_Q )
                    {
                        S.WriteInitialVertexCoordiantes(p_out,k);
                    }
                    
                    if constexpr ( x_out_Q > 0 )
                    {
                        S.WriteInitialEdgeVectors(x_out,k);
                    }
                    
                    if constexpr ( w_Q )
                    {
                        S.WriteShiftVector(w,k);
                    }
                    
                    if constexpr ( y_Q )
                    {
                        S.WriteEdgeVectors(y,k);
                    }
                    
                    if constexpr ( q_Q )
                    {
                        S.WriteVertexPositions(q,k);
                    }
                    
                    if constexpr ( edge_space_Q )
                    {
                        K_edge_space[k] = S.EdgeSpaceSamplingWeight();
                    }
                    
                    if constexpr ( quot_space_Q )
                    {
                        K_quot_space[k] = S.EdgeQuotientSpaceSamplingWeight();
                    }
                }
                
                Time stop = Clock::now();
                
                logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Tools::Duration(start, stop) ) + "." );
                
            },
            thread_count
        );
    }
