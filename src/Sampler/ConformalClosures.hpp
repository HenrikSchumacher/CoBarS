public:

    virtual void ComputeConformalClosures(
        const Real * restrict const x,
              Real * restrict const w,
              Real * restrict const y,
              Real * restrict const K_edge_space,
              Real * restrict const K_quot_space,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        ptic(ClassName()+"::ComputeConformalClosures");
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                Time start = Clock::now();
                
                const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );

                Sampler S ( EdgeLengths().data(), Rho().data(), EdgeCount(), Settings() );

                for( Int k = k_begin; k < k_end; ++k )
                {
                    S.ReadInitialEdgeVectors(x,k);
                    
                    S.ComputeConformalClosure();
                    
                    S.WriteShiftVector(w,k);
                    
                    S.WriteEdgeVectors(y,k);

                    K_edge_space[k] = S.EdgeSpaceSamplingWeight();

                    K_quot_space[k] = S.EdgeQuotientSpaceSamplingWeight();
                }
                
                Time stop = Clock::now();
                
                logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Tools::Duration(start, stop) ) + "." );
                
            },
            thread_count
        );
        
        ptoc(ClassName()+"::ComputeConformalClosures");
    }
