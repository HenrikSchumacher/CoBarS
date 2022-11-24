#pragma once

namespace CycleSampler
{

    template<typename Real, typename Int>
    class RandomVariableBase;
    
    template<int AmbDim, typename Real, typename Int>
    class RandomVariable;
    
    template<typename Real = double, typename Int = long long>
    struct SamplerSettings
    {
        Real tolerance            = std::sqrt(std::numeric_limits<Real>::epsilon());
        Real give_up_tolerance    = 100 * std::numeric_limits<Real>::epsilon();
        Real regularization       = 1.;
        Int  max_iter             = 1000;
        
        Real Armijo_slope_factor  = 0.01;
        Real Armijo_shrink_factor = 0.5;
        Int  max_backtrackings    = 20;
        
        bool use_linesearch       = true;
        
        SamplerSettings() {}
        
        ~SamplerSettings() = default;
        
        SamplerSettings( const SamplerSettings & other )
        :   tolerance(other.tolerance)
        ,   give_up_tolerance(other.give_up_tolerance)
        ,   regularization(other.regularization)
        ,   max_iter(other.max_iter)
        ,   Armijo_slope_factor(other.Armijo_slope_factor)
        ,   Armijo_shrink_factor(other.Armijo_shrink_factor)
        ,   max_backtrackings(other.max_backtrackings)
        ,   use_linesearch(other.use_linesearch)
        {}
        
        void PrintStats() const
        {
            valprint( "tolerance           ", tolerance           , 16 );
            valprint( "give_up_tolerance   ", give_up_tolerance   , 16 );
            valprint( "regularization      ", regularization      , 16 );
            valprint( "max_iter            ", max_iter            , 16 );
            valprint( "Armijo_slope_factor ", Armijo_slope_factor , 16 );
            valprint( "Armijo_shrink_factor", Armijo_shrink_factor, 16 );
            valprint( "max_backtrackings   ", max_backtrackings   , 16 );
            valprint( "use_linesearch      ", use_linesearch      , 16 );
        }
    };
    
    template<
        typename Real    = double,
        typename Int     = long long
    >
    class SamplerBase
    {
        ASSERT_INT(Int);
        ASSERT_FLOAT(Real);
    
    public:
        
        
        using SpherePoints_T = Tensor2<Real,Int>;
        using SpacePoints_T  = Tensor2<Real,Int>;
        using Weights_T      = Tensor1<Real,Int>;
        
        using Setting_T = SamplerSettings<Real,Int>;
        
    protected:
        
        const Int edge_count = 0;
        
        mutable std::mt19937_64 random_engine;
        
        mutable std::normal_distribution<Real> normal_dist {static_cast<Real>(0),static_cast<Real>(1)};
        
        Setting_T settings;
        
        
    public:
        
        SamplerBase()
        {
            std::random_device r;
            
            std::seed_seq seed { r(), r(), r(), r() };
            
//            std::seed_seq seed {1};
            
            random_engine = std::mt19937_64( seed );
        }
        
        explicit SamplerBase(
            const Int edge_count_,
            const Setting_T settings_ = Setting_T()
        )
        :   edge_count(edge_count_)
        ,   settings(settings_)
        {
            std::random_device r;
            
            std::seed_seq seed { r(), r(), r(), r() };
            
//            std::seed_seq seed {1};
            
            random_engine = std::mt19937_64( seed );
        }
        
        virtual ~SamplerBase() = default;
        
    public:
        
        Int EdgeCount() const
        {
            return edge_count;
        }
        
        Int VertexCount() const
        {
            return edge_count + 1;
        }
        
        virtual void Optimize() = 0;
        
        virtual void OptimizeBatch(
                  Real * const x_in,
                  Real * const w_out,
                  Real * const y_out,
            const Int sample_count,
            const Int thread_count = 1,
            bool normalize = true
        ) = 0;
        
//        virtual void SetInitialData(
//            const Real * restrict const x_in,
//            const Real * restrict const omega_in,
//            const bool normalize = true
//        ) = 0;
        
        virtual const Real * InitialEdgeCoordinates() const = 0;
        
        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, bool normalize = true ) = 0;
        
        virtual void ReadInitialEdgeCoordinates( const Real * const x_in, const Int k, bool normalize = true ) = 0;
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out ) const = 0;
        
        virtual void WriteInitialEdgeCoordinates( Real * x_out, const Int k ) const = 0;
        
        
        virtual const Real * EdgeLengths() const = 0;
        
        virtual void ReadEdgeLengths( const Real * const r_in ) = 0;
        
        
        virtual const Real * Rho() const = 0;
        
        virtual void ReadRho( const Real * const rho_in ) = 0;
        
                
        virtual const Real * EdgeCoordinates() const = 0;
        
        virtual void ReadEdgeCoordinates( const Real * const y_in ) = 0;
        
        virtual void ReadEdgeCoordinates( const Real * const y_in, const Int k ) = 0;
        
        virtual void WriteEdgeCoordinates( Real * y_out ) const = 0;
        
        virtual void WriteEdgeCoordinates( Real * y_out, const Int i ) const = 0;
        
        
//        virtual Real * SpaceCoordinates() = 0;
        
        virtual const Real * SpaceCoordinates() const = 0;
        
        virtual void WriteSpaceCoordinates( Real * y_in ) const = 0;
        
        virtual void WriteSpaceCoordinates( Real * y_in, const Int i ) const = 0;
        
        virtual void ComputeSpaceCoordinates() const = 0;
        
        
        
        virtual void ComputeShiftVector() = 0;
        
        virtual void ReadShiftVector( const Real * const w_in ) = 0;
        
        virtual void ReadShiftVector( const Real * const w_in, const Int i ) = 0;
        
        virtual void WriteShiftVector( Real * w_out ) const = 0;
        
        virtual void WriteShiftVector( Real * w_out, const Int i ) const = 0;
        
        virtual Real Residual() const = 0;
        
        virtual Real ErrorEstimator() const = 0;
        
        virtual Int IterationCount() const = 0;
        
        virtual Int MaxIterationCount() const = 0;

        virtual Int AmbientDimension() const = 0;
        
        virtual void RandomSphericalPoints( Real * restrict x_out, const Int m, const Int thread_count = 1 ) const = 0;
        
        virtual void RandomizeInitialEdgeCoordinates() = 0;
        
        
        virtual Real EdgeSpaceSamplingWeight() const = 0;
        
        virtual Real EdgeQuotientSpaceSamplingCorrection() const = 0;
        
        virtual Real EdgeQuotientSpaceSamplingWeight() const = 0;

        virtual void RandomClosedPolygons(
                  Real * const restrict x_out,
                  Real * const restrict w_out,
                  Real * const restrict y_out,
                  Real * const restrict K_edge_space,
                  Real * const restrict K_edge_quotient_space,
            const Int sample_count,
            const Int thread_count = 1
        ) const = 0;
        
        virtual void Sample_Binned(
            Real * restrict bins_out,
            const Int bin_count,
            Real * restrict moments_out,
            const Int moment_count,
            const Real * restrict ranges_out,
            const std::vector< std::unique_ptr<RandomVariableBase<Real,Int>> > & F_list,
            const Int sample_count,
            const Int thread_count = 1
        ) const = 0;
        
        virtual void NormalizeBinnedSamples(
            Real * restrict bins,
            const Int bin_count,
            Real * restrict moments,
            const Int moment_count,
            const Int fun_count
        ) const = 0;
        
        
        const Setting_T & Settings() const
        {
            return settings;
        }
        
    public:
        
        virtual std::string ClassName() const
        {
            return "SamplerBase<"+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
    };
        
} // namespace CycleSampler
