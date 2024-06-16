public:
    
    virtual void ComputeConformalClosure() override
    {
        computeConformalClosure<true,true>();
    }
    
private:
    
    template<bool vertex_pos_Q, bool quot_space_Q>
    void computeConformalClosure()
    {
        ComputeInitialShiftVector();

        Optimize();
        
        if constexpr ( vertex_pos_Q )
        {
            ComputeVertexPositions();
        }
        
        ComputeEdgeSpaceSamplingWeight();
        
        if constexpr ( quot_space_Q )
        {
            ComputeEdgeQuotientSpaceSamplingWeight();
        }
    }
