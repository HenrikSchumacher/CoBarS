#pragma once

namespace CycleSampler
{
    template<int AmbDim, typename Real = double, typename Int = int_fast32_t>
    class SamplerBase
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        using Vector_T          = Tiny::Vector           <AmbDim,Real,Int>;
        using SquareMatrix_T    = Tiny::Matrix           <AmbDim,AmbDim,Real,Int>;
        using SymmetricMatrix_T = Tiny::SelfAdjointMatrix<AmbDim,Real,Int>;
        
        using RandomVariable_T  = RandomVariable<SamplerBase<AmbDim,Real,Int>>;
        
        using Weights_T         = Tensor1<Real,Int>;
        using Setting_T         = SamplerSettings<Real,Int>;
        
    protected:
        
        Setting_T settings;
        
    public:
        
        SamplerBase( SamplerSettings<Real,Int> settings_ )
        :   settings(  settings_ )
        {}
        
        virtual ~SamplerBase() = default;
        
        virtual void Optimize() = 0;
        
        virtual Int EdgeCount() const = 0;
        
        virtual Real EdgeSpaceSamplingWeight() const = 0;
        
        virtual Real EdgeQuotientSpaceSamplingCorrection() const = 0;
        
        virtual Real EdgeQuotientSpaceSamplingWeight() const = 0;
        
//        virtual const SpherePoints_T & InitialEdgeCoordinates() const = 0;
        
        virtual Real InitialEdgeCoordinates( const Int k, const Int i ) const = 0;
        
        virtual Vector_T InitialEdgeCoordinates( const Int k ) const = 0;
        
        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, bool normalize = true ) = 0;
        
        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, const Int k, bool normalize = true ) = 0;
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out ) const = 0;
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out, const Int k ) const = 0;
        
        virtual void RandomizeInitialEdgeCoordinates() = 0;
        
//        virtual const SpherePoints_T & EdgeCoordinates() const = 0;
        
        virtual Real EdgeCoordinates( const Int k, const Int i ) const = 0;
        
        virtual Vector_T EdgeCoordinates( const Int k ) const = 0;
        
        virtual void ReadEdgeCoordinates( ptr<Real> y_in ) = 0;
        
        virtual void ReadEdgeCoordinates( ptr<Real> y_in, const Int k ) = 0;
        
        virtual void WriteEdgeCoordinates( mut<Real> y_out ) const = 0;
        
        virtual void WriteEdgeCoordinates( Real * y_out, const Int k ) const = 0;
        
//        virtual const SpacePoints_T & SpaceCoordinates() const = 0;
        
        virtual Real SpaceCoordinates( const Int k, const Int i ) const = 0;
        
        virtual Vector_T SpaceCoordinates( const Int k ) const = 0;
        
        virtual void WriteSpaceCoordinates( Real * p_out ) const = 0;
        
        virtual void WriteSpaceCoordinates( Real * p_out, const Int k ) const = 0;
        
        virtual void ComputeSpaceCoordinates() = 0;
        
        
        virtual const Weights_T & EdgeLengths() const = 0;
        
        virtual void ReadEdgeLengths( const Real * const r_in ) = 0;
        
        virtual const Weights_T & Rho() const = 0;
        
        virtual void ReadRho( const Real * const rho_in ) = 0;
        
        virtual void ComputeShiftVector() = 0;
        
        virtual  void ReadShiftVector( const Real * const w_in ) = 0;
        
        virtual void ReadShiftVector( const Real * const w_in, const Int k ) = 0;
        
        virtual void WriteShiftVector( Real * w_out ) const = 0;
        
        virtual void WriteShiftVector( Real * w_out, const Int k ) const = 0;
        
        virtual const Vector_T & ShiftVector() const = 0;
        
        virtual Real Residual() const = 0;
        
        virtual Real ErrorEstimator() const = 0;
        
        virtual Int IterationCount() const = 0;
        
        virtual Int MaxIterationCount() const = 0;
        
        virtual void OptimizeBatch(
                                   ptr<Real>  x_in,
                                   mut<Real>  w_out,
                                   mut<Real>  y_out,
                                   const Int  sample_count,
                                   const Int  thread_count = 1,
                                   const bool normalize = true
                                   ) = 0;
        
        virtual void RandomClosedPolygons(
                                          mut<Real> x_out,
                                          mut<Real> w_out,
                                          mut<Real> y_out,
                                          mut<Real> K_edge_space,
                                          mut<Real> K_edge_quotient_space,
                                          const Int sample_count,
                                          const Int thread_count = 1
                                          ) const = 0;
        
        virtual void Sample(
                            mut<Real> sampled_values,
                            mut<Real> edge_space_sampling_weights,
                            mut<Real> edge_quotient_space_sampling_weights,
                            std::shared_ptr<RandomVariable_T> & F_,
                            const Int sample_count,
                            const Int thread_count = 1
                            ) const = 0;
        
        virtual void Sample(
                            mut<Real> sampled_values,
                            mut<Real> edge_space_sampling_weights,
                            mut<Real> edge_quotient_space_sampling_weights,
                            const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
                            const Int sample_count,
                            const Int  thread_count = 1
                            ) const = 0;
        
        virtual void SampleCompressed(
                                      mut<Real> bins_out,
                                      const Int bin_count_,
                                      mut<Real> moments_out,
                                      const Int moment_count_,
                                      ptr<Real> ranges,
                                      const std::vector< std::shared_ptr<RandomVariable_T> > & F_list_,
                                      const Int sample_count,
                                      const Int thread_count = 1
                                      ) const = 0;
        
        const Setting_T & Settings() const
        {
            return settings;
        }
        
        Int AmbientDimension() const
        {
            return AmbDim;
        }
        
        
        void NormalizeCompressedSamples(
                                        mut<Real> bins,
                                        const Int bin_count,
                                        mut<Real> moments,
                                        const Int moment_count,
                                        const Int fun_count
                                        ) const
        {
            ptic(this->ClassName()+"::NormalizeCompressedSamples");
            for( Int i = 0; i < 3; ++i )
            {
                for( Int j = 0; j < fun_count; ++j )
                {
                    // Normalize bins and moments.
                    
                    mut<Real> bins_i_j = &bins[ (i*fun_count+j)*bin_count ];
                    
                    mut<Real> moments_i_j = &moments[ (i*fun_count+j)*moment_count ];
                    
                    // The field for zeroth moment is assumed to contain the total mass.
                    Real factor = Real(1)/moments_i_j[0];
                    
                    scale_buffer( factor, bins_i_j,    bin_count    );
                    
                    scale_buffer( factor, moments_i_j, moment_count );
                }
            }
            ptoc(this->ClassName()+"::NormalizeCompressedSamples");
        }
        
        virtual std::string ClassName() const
        {
            return std::string("SamplerBase") + "<" + ToString(AmbDim) + "," + TypeName<Real> + "," + TypeName<Int>  + ">";
        }
    };
    
} // namespace CycleSampler
