public:

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
            CreatePolygons<0,0,0,0,0,1,0,0,1>(
                nullptr, nullptr, nullptr, nullptr, nullptr, y, nullptr, nullptr, K,
                sample_count, thread_count
            );
        }
        else
        {
            CreatePolygons<0,0,0,0,0,1,0,1,0>(
                nullptr, nullptr, nullptr, nullptr, nullptr, y, nullptr, K, nullptr,
                sample_count, thread_count
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

        CreatePolygons<0,1,0,0,1,1,0,1,1>(
            nullptr, x, nullptr, nullptr, w, y, nullptr, K_edge_space, K_quot_space,
            sample_count, thread_count
        );
        
        ptoc(ClassName()+"::CreateRandomCentralizedPointClouds_Detailed");
    }

