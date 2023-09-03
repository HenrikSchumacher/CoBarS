protected:
    
    const Int edge_count = 0;
    
    mutable PRNG_T random_engine [AmbDim];
    
    mutable std::normal_distribution<Real> normal_dist {zero,one};
    
    Setting_T settings;

    std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> x;
    std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> y;
    std::conditional_t<vectorizeQ, VectorList_T, Matrix_T> p;
    
    Weights_T      r {0};
    Weights_T    rho {0};
    
    Real total_r_inv = one;
    
    Vector_T w;           // current point in hyperbolic space.
    Vector_T F;           // right hand side of Newton iteration.
    SymmetricMatrix_T DF; // nabla F(0) with respect to measure ys
    SymmetricMatrix_T L;  // storing Cholesky factor.
    Vector_T u;           // update direction
    
    Vector_T z;           // Multiple purpose buffer.
    
    Int iter = 0;
    
    Real squared_residual = 1;
    Real         residual = 1;
    
    Real edge_space_sampling_weight = 0;
    Real edge_quotient_space_sampling_correction = 0;
    
    Real lambda_min = eps;
    Real q = one;
    Real errorestimator = infty;
    
    bool linesearchQ = true;    // Toggle line search.
    bool succeededQ  = false;   // Whether algorithm has succeded.
    bool continueQ   = true;    // Whether to continue with the main loop.
    bool ArmijoQ     = false;   // Whether Armijo condition was met last time we checked.
    
    static constexpr Real zero              = 0;
    static constexpr Real half              = 0.5;
    static constexpr Real one               = 1;
    static constexpr Real two               = 2;
    static constexpr Real three             = 3;
    static constexpr Real four              = 4;
    static constexpr Real eps               = std::numeric_limits<Real>::eps();
    static constexpr Real infty             = std::numeric_limits<Real>::max();
    static constexpr Real small_one         = 1 - 16 * eps;
    static constexpr Real big_one           = 1 + 16 * eps;
    static constexpr Real g_factor          = 4;
    static constexpr Real g_factor_inv      = one/g_factor;
    static constexpr Real norm_threshold    = static_cast<Real>(0.99 * 0.99 + 16 * eps);
    static constexpr Real two_pi            = Scalar::TwoPi<Real>;
