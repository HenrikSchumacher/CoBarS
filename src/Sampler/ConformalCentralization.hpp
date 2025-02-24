public:

    virtual void ComputeConformalCentralizations(
        const Real * restrict const x,
              Real * restrict const w,
              Real * restrict const y,
              Real * restrict const K_edge_space,
              Real * restrict const K_quot_space,
        const Int sample_count,
        const Int thread_count = 1
    ) const override
    {
        TOOLS_PTIC(ClassName()+"::ComputeConformalCentralizations");
        
        CreatePolygons<1,0,0,0,1,1,0,1,1>(
            x, nullptr, nullptr, nullptr, w, y, nullptr, K_edge_space, K_quot_space,
            sample_count, thread_count
        );
        
        TOOLS_PTOC(ClassName()+"::ComputeConformalCentralizations");
    }
