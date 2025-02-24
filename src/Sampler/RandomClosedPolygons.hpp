public:

    virtual void CreateRandomClosedPolygons(
        Real * restrict const q,
        Real * restrict const K,
        const Int sample_count,
        const bool quotient_space_Q = true,
        const Int thread_count = 1
    ) const override
    {
        TOOLS_PTIC(ClassName()+"::CreateRandomClosedPolygons");
        
        if( quotient_space_Q )
        {
            CreatePolygons<0,0,0,0,0,0,1,0,1>(
                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, q, nullptr, K,
                sample_count, thread_count
            );
        }
        else
        {
            CreatePolygons<0,0,0,0,0,0,1,1,0>(
                nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, q, K, nullptr,
                sample_count, thread_count
            );
        }
    
        TOOLS_PTOC(ClassName()+"::CreateRandomClosedPolygons");
    }
