public:

    virtual void ConformalClosures(
        ptr<Real> x_in,
        mut<Real> w_out,
        mut<Real> y_out,
        mut<Real> K_edge_space,
        mut<Real> K_edge_quotient_space,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        ptic(ClassName()+"::ConformalClosures");
        
        ParallelDo(
            [&,this]( const Int thread )
            {
                Time start = Clock::now();
                
                const Int k_begin = JobPointer( sample_count, thread_count, thread     );
                const Int k_end   = JobPointer( sample_count, thread_count, thread + 1 );

                Sampler S ( EdgeLengths().data(), Rho().data(), edge_count, Settings() );

                for( Int k = k_begin; k < k_end; ++k )
                {
                    S.ReadInitialEdgeCoordinates(x_in, k);
                    
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
        
        ptoc(ClassName()+"::ConformalClosures");
    }
