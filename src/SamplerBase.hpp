#pragma once

namespace CoBarS
{
    /*!
     * @brief Serves as base class to provide runtime polymorphism for the template classes CoBarS::Sampler
     *
     * @tparam AmbDim_ The dimension of the ambient space.
     *
     * @tparam Real_ A real floating point type.
     *
     * @tparam Int_  An integer type.
     */
    
    template<
        int AmbDim_,
        typename Real_ = double,
        typename Int_  = int_fast32_t
    >
    class SamplerBase
    {
        static_assert(FloatQ<Real_>,"");
        static_assert(Scalar::RealQ<Real_>,"");
        static_assert(IntQ<Int_>,"");
        
    public:
        
        using Real = Real_;
        using Int  = Int_;
        
        static constexpr Int AmbDim = int_cast<Int>(AmbDim_);
        
        using Vector_T          = Tiny::Vector           <AmbDim,Real,Int>;
        using SquareMatrix_T    = Tiny::Matrix           <AmbDim,AmbDim,Real,Int>;
        using SymmetricMatrix_T = Tiny::SelfAdjointMatrix<AmbDim,Real,Int>;
        
        using RandomVariable_T  = RandomVariable<SamplerBase<AmbDim,Real,Int>>;
        
        using Weights_T         = Tensor1<Real,Int>;
        using Setting_T         = SamplerSettings<Real,Int>;
        
    protected:
        
        Setting_T settings_;
        
    public:
        
        /*!
         * @brief Use this constructor if you want to use user-defined options. 
         */
        
        SamplerBase( SamplerSettings<Real,Int> settings )
        :   settings_(  settings )
        {}
        
        SamplerBase() = default;
        
        virtual ~SamplerBase() = default;
        
    public:
        
        /*!
         * @brief Returns the number of edges.
         */
        
        virtual Int EdgeCount() const = 0;
        
        /*!
         * @brief After the conformal closure has been computed, this returns the reweighting factor to be used for sampling on the space of closed polygons _without_ modding-out the rotation group.
         */
        
        virtual Real EdgeSpaceSamplingWeight() const = 0;
        
        /*!
         * @brief After the conformal closure has been computed, this returns the reweighting factor to be used for sampling on the space of closed polygons modulo the rotation group.
         */
        
        virtual Real EdgeQuotientSpaceSamplingWeight() const = 0;
        
//        /*!
//         * @brief Returns the `i`-th coordinate of the `k`-th edge vector of the open polygon. This routine is slow; try to avoid it.
//         */
//        
//        virtual Real InitialEdgeCoordinate( const Int k, const Int i ) const = 0;
        
        /*!
         * @brief Returns the the `k`-th edge vector of the open polygon.
         */
        
        virtual Vector_T InitialEdgeVector( const Int k ) const = 0;
        
        /*!
         * @brief Reads a set of edge vectors for the closed polygon from the buffer `x`.
         *
         * @param x Source array; assumed to be of size of at least `EdgeCount() * AmbientDimension() * (offset + 1)`. It is assumed that that the coordinates are stored in interleaved form and that the `EdgeCount() * AmbientDimension()` coordinates of a single polygon are stored consecutively.
         *
         * @param offset Indicates that the polygon starts at `&[EdgeCount() * AmbientDimension() * offset]`.
         *
         * @param normalizeQ Only if `normalizeQ` is set to true, the vectors get normalized.
         */
        
        virtual void ReadInitialEdgeVectors( const Real * restrict const x, const Int offset = 0, bool normalizeQ = true ) = 0;

        /*!
         * @brief Writes the edge vectors of the open polygon to buffer `x`.
         *
         * @param x Target array; assumed to be of size of at least `EdgeCount() * AmbientDimension() * (offset + 1)`. The coordinates are stored in interleaved form and the `EdgeCount() * AmbientDimension()` coordinates of a single polygon are stored consecutively.
         *
         * @param offset Indicates that the polygon starts at `&[EdgeCount() * AmbientDimension() * offset]`.
         */
        
        virtual void WriteInitialEdgeVectors( Real * restrict const x, const Int offset = 0 ) const = 0;
        
        /*!
         * @brief Fills the internal buffer of the open polygon with random unit vectors, effectively generating a random open polygon.
         */
        
        virtual void RandomizeInitialEdgeVectors() = 0;
        
//        /*!
//         * @brief Returns the `i`-th coordinate of the `k`-th edge vector of the _closed_ polygon. This routine is slow; try to avoid it.
//         */
//        
//        virtual Real EdgeCoordinate( const Int k, const Int i ) const = 0;
        
        /*!
         * @brief Returns the coordinate of the `k`-th edge vector of the _closed_ polygon.
         */
        
        virtual Vector_T EdgeVector( const Int k ) const = 0;
        
        /*!
         * @brief Writes the edge vectors of the closed polygon to buffer `y`.
         *
         * @param y Target array; assumed to be of size of at least `EdgeCount() * AmbientDimension() * (offset + 1)`. The coordinates are stored in interleaved form and the `EdgeCount() * AmbientDimension()` coordinates of a single polygon are stored consecutively.
         *
         * @param offset Indicates that the polygon starts at `&[EdgeCount() * AmbientDimension() * offset]`.
         */
        
        virtual void WriteEdgeVectors( Real * restrict const y, const Int offset = 0 ) const = 0;
        

        /*!
         * @brief Returns the position of `k`-th edge vertex of the closed polygon.
         */
        
        virtual Vector_T VertexPosition( const Int k ) const = 0;
        
//        virtual Real VertexCoordinate( const Int k, const Int i ) const = 0;
        
        /*!
         * @brief Writes the vertex positions of the closed polygon to buffer `p`.
         *
         * @param p Target array; assumed to be of size of at least `(EdgeCount()+1) * AmbientDimension() * (offset + 1)`. The coordinates are stored in interleaved form and the `EdgeCount() * AmbientDimension()` coordinates of a single polygon are stored consecutively. IMPORTANT: The position of the first vertex is appended to the list, so this are `EdgeCount()+1` vertices.
         *
         * @param offset Indicates that the polygon starts at `&[EdgeCount() * AmbientDimension() * offset]`.
         */
        
        virtual void WriteVertexPositions( Real * restrict const p, const Int offset = 0 ) const = 0;
        
        /*!
         * @brief Makes sure that the vertex coordinates of the closed polygon are computed.
         */
        
    private:
        
        virtual void ComputeVertexPositions() const = 0;

        /*!
         * @brief Returns the list of edge lengths.
         */
        
    public:
        
        virtual const Weights_T & EdgeLengths() const = 0;
        
        /*!
         * @brief Reads a new list of edge lengths from buffer `r`.
         *
         * @param r Buffer containing the new edge lengths; assumed to have length `EdgeCount()`.
         */
        
        virtual void ReadEdgeLengths( const Real * const r ) = 0;
        
        /*!
         * @brief Returns the list of weights for the Riemannian metrics on the product of unit spheres.
         */
        
        virtual const Weights_T & Rho() const = 0;
        
        /*!
         * @brief Reads a new list of edge weights for the Riemannian metric from buffer `rho`.
         *
         * @param rho Buffer containing the new weights; assumed to have length `EdgeCount()`.
         */
        
        virtual void ReadRho( const Real * const rho ) = 0;
        
        /*!
         * @brief Compute an initial guess for the conformal barycenter to be used in the optimization routine `Optimize`.
         */
        
        virtual void ComputeInitialShiftVector() = 0;
        
        /*!
         * @brief Reads an initial guess for the conformal barycenter (to be used in the optimization routine `Optimize`) from the position `AmbientDimension() * offset` in the buffer `w`.
         */
        
        virtual void ReadShiftVector( const Real * const w, const Int offset = 0 ) = 0;
        

        
        /*!
         * @brief Writes the current shift vector (e.g., the conformal barycenter after the optimization has succeeded) to the position `AmbientDimension() * offset` in the buffer `w`.
         */
        
        virtual void WriteShiftVector( Real * w, const Int offset = 0 ) const = 0;
        
        /*!
         * @brief Returns the current shift vector (e.g., the conformal barycenter after the optimization has succeeded).
         */
        
        virtual const Vector_T & ShiftVector() const = 0;
        
        /*!
         * @brief Returns the current residual of the optimization routine `Optimize`.
         */
        
        virtual Real Residual() const = 0;
        
        /*!
         * @brief After `Optimize` succeeded, this is an upper bound for the distance of `ShiftVector()` to the true conformal barycenter of the open polygon.
         */
        
        virtual Real ErrorEstimator() const = 0;
        
        /*!
         * @brief Returns the number of iterations the last call to `Optimize` needed.
         */
        
        virtual Int IterationCount() const = 0;
        
        /*!
         * @brief Returns the maximal number of iterations `Optimize` is allowed to take.
         */
        
        virtual Int MaxIterationCount() const = 0;
        
        /*!
         * @brief Uses the conformal closure procedure to close the open polygon (loaded with `ReadInitialEdgeVectors` or generated with `ReadInitialEdgeVectors`). Afterwards, the resulting closed polygon can be accessed with routines `*EdgeVectors`. The conformal barycenter can be accessed with `*ShiftVector`.
         */
        
        virtual void ComputeConformalClosure() = 0;
        

        /*!
         * @brief Generates `sample_count` random open polygons and writes them to the buffer `x`.
         *
         * @param x The output array for the open polygons; it is assumed to have size at least `sample_count * EdgeCount() * AmbientDimension()`.
         *
         * @param sample_count Number of polygons to generated.
         *
         * @param thread_count Number of threads to use.
         */
        
        virtual void RandomOpenPolygons(
            Real * restrict const x,
            const Int sample_count,
            const Int thread_count = 1
        ) const = 0;
        
        /*!
         * @brief Generates `sample_count` random open polygons, closes them, and writes the relevant information to the supplied buffers.
         *
         * @param x The output array for the _open_ polygons; it is assumed to have size at least `sample_count * EdgeCount() * AmbientDimension()`.
         *
         * @param w The output array for the conformal barycenters polygons; it is assumed to have size at least `sample_count * AmbientDimension()`.
         *
         * @param y The output array for the _closed_ polygons; it is assumed to have size at least `sample_count * EdgeCount() * AmbientDimension()`.
         *
         * @param K_edge_space The output array for the reweighting factors of the polygons space (rotation group not modded out); it is assumed to have size at least `sample_count`.
         *
         * @param K_quot_space The output array for the reweighting factors of the polygons space modulo rotation group; it is assumed to have size at least `sample_count`.
         *
         * @param sample_count Number of polygons to generated.
         *
         * @param thread_count Number of threads to use.
         */
        
        virtual void RandomClosedPolygons(
            Real * restrict const x,
            Real * restrict const w,
            Real * restrict const y,
            Real * restrict const K_edge_space,
            Real * restrict const K_quot_space,
            const Int sample_count,
            const Int thread_count = 1
        ) const = 0;
        
        /*!
         * @brief Reads `sample_count` open polygons from buffer `x`, closes them, and writes the relevant information to the supplied buffers.
         *
         * @param x The input array for the _open_ polygons; it is assumed to have size at least `sample_count * EdgeCount() * AmbientDimension()`.
         *
         * @param w The output array for the conformal barycenters polygons; it is assumed to have size at least `sample_count * AmbientDimension()`.
         *
         * @param y The output array for the _closed_ polygons; it is assumed to have size at least `sample_count * EdgeCount() * AmbientDimension()`.
         *
         * @param K_edge_space The output array for the reweighting factors of the polygons space (rotation group not modded out); it is assumed to have size at least `sample_count`.
         *
         * @param K_quot_space The output array for the reweighting factors of the polygons space modulo rotation group; it is assumed to have size at least `sample_count`.
         *
         * @param sample_count Number of polygons to generated.
         *
         * @param thread_count Number of threads to use.
         */
        
        virtual void ConformalClosures(
            const Real * restrict const x,
                  Real * restrict const w,
                  Real * restrict const y,
                  Real * restrict const K_edge_space,
                  Real * restrict const K_quot_space,
            const Int sample_count,
            const Int thread_count = 1
        ) const = 0;
        
        
        /*!
         * @brief Generates `sample_count` random open polygons, closes them, evaluated the random variable `F` on them, and writes the relevant information to the supplied buffers. The generated polygons are discarded immediately after evaluating the random variables on them.
         *
         * @param sampled_values The output array for sampled function values; expected to be an array of size at least `sample_count`.
         *
         * @param K_edge_space The output array for the reweighting factors of the polygons space (rotation group not modded out); it is assumed to have size at least `sample_count`.
         *
         * @param K_quot_space The output array for the reweighting factors of the polygons space modulo rotation group; it is assumed to have size at least `sample_count`.
         * 
         * @param F The random variable to sample.
         *
         * @param sample_count Number of polygons to generated.
         *
         * @param thread_count Number of threads to use.
         */
        
        virtual void Sample(
            Real * restrict const sampled_values,
            Real * restrict const K_edge_space,
            Real * restrict const K_quot_space,
            std::shared_ptr<RandomVariable_T> & F,
            const Int sample_count,
            const Int thread_count = 1
        ) const = 0;
        
        /*!
         * @brief Generates `sample_count` random open polygons, closes them, evaluated some random variables on them, and writes the relevant information to the supplied buffers. The generated polygons are discarded immediately after evaluating the random variables on them.
         *
         * @param sampled_values The output array for sampled function values; expected to be an array of size at least `sample_count * fun_count`. Samples are stored in interleaved form, the function values for each polygon are contiguous in this buffer.
         *
         * @param K_edge_space The output array for the reweighting factors of the polygons space (rotation group not modded out); it is assumed to have size at least `sample_count`.
         *
         * @param K_quot_space The output array for the reweighting factors of the polygons space modulo rotation group; it is assumed to have size at least `sample_count`.
         *
         * @param F_list The random variables to sample.
         *
         * @param sample_count Number of polygons to generated.
         *
         * @param thread_count Number of threads to use.
         */
        
        virtual void Sample(
            Real * restrict const sampled_values,
            Real * restrict const K_edge_space,
            Real * restrict const K_quot_space,
            const std::vector< std::shared_ptr<RandomVariable_T> > & F_list,
            const Int sample_count,
            const Int thread_count = 1
        ) const = 0;
        
        
        /*!
         * @brief Returns the list of loaded random variables.
         */
        
        virtual const std::vector<std::shared_ptr<RandomVariable_T>> & RandomVariables() const = 0;
        
        /*!
         * @brief Returns the number of loaded random variables.
         */
        
        virtual Int RandomVariablesCount() const = 0;
        
        /*!
         * @brief Loads a list of random variables.
         */
        
        virtual void LoadRandomVariables( const std::vector<std::shared_ptr<RandomVariable_T>> & F_list ) const = 0;
        
        /*!
         * @brief Clears the list of random variables.
         */
        
        virtual void ClearRandomVariables() const = 0;
        
        /*!
         * @brief Returns the `i`-th loaded random variable, evaluated on the current `EdgeVectors`.
         */
        
        virtual Real EvaluateRandomVariable( Int i ) const = 0;

        
        /*!
         * @brief Returns the current setting of the `CoBarS::SamplerBase` instance.
         */
        
        const Setting_T & Settings() const
        {
            return settings_;
        }
        
        /*!
         * @brief Returns the dimension of the ambient space.
         */
        
        Int AmbientDimension() const
        {
            return AmbDim;
        }
        
        
        /*!
         * @brief Generates `sample_count` random closed polygons, evaluates the random variables, and writes the results to the supplied buffers.
         *
         * @param bins Buffer that contains the bins. Assumed to be of size `3 * fun_count * bin_count`. There are three slice of dimension `fun_count * bin_count`, each representing the bins for (0) unreweighted samples, (1) samples reweighted w.r.t to probability density of polygon space (without modding out rotation group), and (2) samples reweighted w.r.t to probability density of polygon spaces modulo rotation group, respectively.
         *
         * @param bin_count The number of bins used for binning.
         *
         * @param moms Buffer for the moments. Assumed to be of size `3 * fun_count * mom_count`. There are three slice of dimension `fun_count * mom_count`, each representing the moments for (0) unreweighted samples, (1) samples reweighted w.r.t to probability density of polygon space (without modding out rotation group), and (2) samples reweighted w.r.t to probability density of polygon spaces modulo rotation group, respectively.
         *
         * @param mom_count The number of moments to compute.
         *
         * @param random_vars The list of random variables to sample.
         */
        
        virtual void BinnedSample(
                  Real * restrict const bins,    const Int bin_count,
                  Real * restrict const moms, const Int mom_count,
            const Real * restrict const ranges,
            const std::vector< std::shared_ptr<RandomVariable_T> > & random_vars,
            const Int sample_count,
            const Int thread_count = 1
        ) const = 0;
        
        
        /*!
         * @brief Normalizes samples generated by the `BinnedSample` routine.
         *
         * @param bins Buffer that contains the bins. Assumed to be of size `3 * fun_count * bin_count`. There are three slice of dimension `fun_count * bin_count`, each representing the bins for (0) unreweighted samples, (1) samples reweighted w.r.t to probability density of polygon space (without modding out rotation group), and (2) samples reweighted w.r.t to probability density of polygon spaces modulo rotation group, respectively.
         *
         * @param bin_count The number of bins used for binning.
         *
         * @param moms Buffer for the moments. Assumed to be of size `3 * fun_count * mom_count`. There are three slice of dimension `fun_count * mom_count`, each representing the moments for (0) unreweighted samples, (1) samples reweighted w.r.t to probability density of polygon space (without modding out rotation group), and (2) samples reweighted w.r.t to probability density of polygon spaces modulo rotation group, respectively.
         *
         * @param mom_count The number of moments to compute.
         *
         * @param random_vars The list of random variables to sample.
         */
        
        void NormalizeBinnedSamples(
            Real * restrict const  bins, const Int bin_count,
            Real * restrict const  moms, const Int mom_count,
            const std::vector< std::shared_ptr<RandomVariable_T> > & random_vars
        ) const
        {
            const Int f_count = static_cast<Int>(random_vars.size());
            
            ptic(this->ClassName()+"::BinnedSamples");
            for( Int i = 0; i < 3; ++i )
            {
                for( Int j = 0; j < f_count; ++j )
                {
                    // Normalize bins and moments.
                    
                    mptr<Real> bins_i_j = &bins[ (i*f_count+j)*bin_count ];
                    
                    mptr<Real> moms_i_j = &moms[ (i*f_count+j)*mom_count ];
                    
                    // The field for zeroth moment is assumed to contain the total mass.
                    Real factor = Real(1)/moms_i_j[0];
                    
                    scale_buffer( factor, bins_i_j, bin_count    );
                    
                    scale_buffer( factor, moms_i_j, mom_count );
                }
            }
            ptoc(this->ClassName()+"::BinnedSamples");
        }
        
        
        /*!
         * @brief Generates random closed polygons and samples the mean of random variables on them, until the radii of the confidence interval to level `confidence_level` are small enough.
         *
         * @param random_vars The list of random variables to sample.
         *
         * @param means A buffer for storing the means; assumed to be of size `random_vars.size()`.
         *
         * @param sample_variances A buffer for storing the sampling variances; assumed to be of size `random_vars.size()`.
         *
         * @param errors A buffer for storing the sampling errors; assumed to be of size `random_vars.size()`.
         *
         * @param radii A buffer for storing the tolerances (desired radii of the confidence intervals); assumed to be of size `random_vars.size()`.
         *
         * @param max_sample_count Maximum number of samples to take.
         *
         * @param quotient_space_Q Whether to do the sampling in the quotient space w.r.t. to the action of the rotation group (`quotient_space_Q == true`) or not (`quotient_space_Q == false`).
         *
         * @param thread_count Number of threads to use.
         *
         * @param confidence_level The confidence level to use to compute the confidence intervals.
         *
         * @param chunk_size How many random polygons ought to be sampled at once per thread. Larger values need more memory. Slower values reduce the overhead a bit.
         *
         * @param relativeQ Whether the `radii` are to be interpreted as relative to the (unknown) mean (`relativeQ == true`) or not (`relativeQ == false`).
         *
         * @param verboseQ Whether to print some intermediate information (`verboseQ == true`) or not (`verboseQ == false`).
         *
         * @return The total number of samples needed.
         */
        
        virtual Int ConfidenceSample(
            const std::vector< std::shared_ptr<RandomVariable_T> > & random_vars,
                  Real * restrict const means,
                  Real * restrict const sample_variances,
                  Real * restrict const errors,
            const Real * restrict const radii,
            const Int  max_sample_count,
            const bool quotient_space_Q,
            const Int  thread_count = 1,
            const Real confidence_level = 0.95,
            const Int  chunk_size = 1000000,
            const bool relativeQ = false,
            const bool verboseQ = true
        ) const = 0;
        
        
        /*! @brief Returns a string that identifies the class' pseudorandom number generator. Good for debugging and printing messages. */
        
        virtual std::string PRNG_Name() const = 0;
        
        /*! @brief Returns a string that identifies the class. Good for debugging and printing messages. */
        
        virtual std::string ClassName() const
        {
            return std::string("CoBarS::SamplerBase") + "<" + ToString(AmbDim) + "," + TypeName<Real> + "," + TypeName<Int>  + ">";
        }
    };
    
} // namespace CoBarS
