public:

    virtual void CreateRandomClosedPolygons(
        Real * restrict const p,
        Real * restrict const K,
        const Int sample_count,
        const bool quotient_space_Q = true,
        const Int thread_count = 1
    ) const override
    {
        ptic(ClassName()+"::CreateRandomClosedPolygons");
        
        if( quotient_space_Q )
        {
            createRandomClosedPolygons<0,0,0,1,0,1>(
                 nullptr, nullptr, nullptr, p, nullptr, K, sample_count, thread_count
            );
        }
        else
        {
            createRandomClosedPolygons<0,0,0,1,1,0>(
                 nullptr, nullptr, nullptr, p, K, nullptr, sample_count, thread_count
            );
        }
    
        ptoc(ClassName()+"::CreateRandomClosedPolygons");
    }


    virtual void CreateRandomCentralizedPointClouds(
        Real * restrict const y,
        Real * restrict const K,
        const Int sample_count,
        const bool quotient_space_Q = true,
        const Int thread_count = 1
    ) const override
    {
        ptic(ClassName()+"::CreateRandomCentralizedPointClouds");

        if( quotient_space_Q )
        {
            createRandomClosedPolygons<0,0,1,0,0,1>(
                 nullptr, nullptr, y, nullptr, nullptr, K, sample_count, thread_count
            );
        }
        else
        {
            createRandomClosedPolygons<0,0,1,0,1,0>(
                 nullptr, nullptr, y, nullptr, K, nullptr, sample_count, thread_count
            );
        }
        
        ptoc(ClassName()+"::CreateRandomCentralizedPointClouds");
    }

    virtual void CreateRandomCentralizedPointClouds_Detailed(
        Real * restrict const x,
        Real * restrict const w,
        Real * restrict const y,
        Real * restrict const K_edge_space,
        Real * restrict const K_quot_space,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        ptic(ClassName()+"::CreateRandomCentralizedPointClouds_Detailed");
  
        createRandomClosedPolygons<1,1,1,0,1,1>(
             x, w, y, nullptr, K_edge_space, K_quot_space, sample_count, thread_count
        );
        
        ptoc(ClassName()+"::CreateRandomCentralizedPointClouds_Detailed");
    }

private:

    // Common implementation of CreateRandomClosedPolygons and CreateRandomCentralizedPointClouds; only to be called from there.
    template<bool x_Q, bool w_Q, bool y_Q, bool p_Q, bool edge_space_Q, bool quot_space_Q>
    void createRandomClosedPolygons(
        mptr<Real> x,
        mptr<Real> w,
        mptr<Real> y,
        mptr<Real> p,
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
                    S.RandomizeInitialEdgeVectors();
                    
                    if constexpr ( x_Q )
                    {
                        S.WriteInitialEdgeVectors(x,k);
                    }
                    
                    S.computeConformalClosure<p_Q,quot_space_Q>();
                    
                    if constexpr ( w_Q )
                    {
                        S.WriteShiftVector(w,k);
                    }
                    
                    if constexpr ( y_Q )
                    {
                        S.WriteEdgeVectors(y,k);
                    }
                    
                    if constexpr ( p_Q )
                    {
                        S.WriteVertexPositions(p,k);
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
