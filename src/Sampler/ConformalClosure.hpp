public:

    virtual void ComputeConformalClosures(
        const Real * restrict const p,
              Real * restrict const w,
              Real * restrict const q,
              Real * restrict const K_edge_space,
              Real * restrict const K_quot_space,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        ptic(ClassName()+"::ComputeConformalClosures");
        
        CreatePolygons<0,0,1,0,1,0,1,1,1>(
            nullptr, nullptr, p, nullptr, w, nullptr, q, K_edge_space, K_quot_space,
            sample_count, thread_count
        );
        
        ptoc(ClassName()+"::ComputeConformalClosures");
    }
