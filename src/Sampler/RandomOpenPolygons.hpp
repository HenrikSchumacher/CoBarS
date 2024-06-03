public:

    virtual void RandomOpenPolygons(
        Real * restrict const x,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        ptic(ClassName()+"::RandomOpenPolygons");
        
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
                    
                    S.WriteInitialEdgeVectors(x,k);
                }
                
                Time stop = Clock::now();
                
                logprint("Thread " + ToString(thread) + " done. Time elapsed = " + ToString( Duration(start, stop) ) + "." );
                
            },
            thread_count
        );
        
        ptoc(ClassName()+"::RandomOpenPolygons");
    }
