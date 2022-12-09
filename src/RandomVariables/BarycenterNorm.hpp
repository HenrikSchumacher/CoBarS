#pragma once

namespace CycleSampler
{
        
#define CLASS BarycenterNorm
        
        template<int AmbDim, typename Real = double, typename Int = long long>
        class CLASS : public RandomVariable<AmbDim,Real,Int>
        {
        private:
            
            using Base_T            = RandomVariable<AmbDim,Real,Int>;
            
        public:
            
            using Sampler_T         = typename Base_T::Sampler_T;
            using SpherePoints_T    = typename Base_T::SpherePoints_T;
            using SpacePoints_T     = typename Base_T::SpacePoints_T;
            using Weights_T         = typename Base_T::Weights_T;
            using Vector_T          = typename Base_T::Vector_T;
            
            CLASS() = default;
            
            virtual ~CLASS() override = default;
            
            __ADD_CLONE_CODE__(CLASS)
            
        protected:
            
            
            virtual Real operator()( const Sampler_T & C ) const override
            {
                const SpacePoints_T & p = C.SpaceCoordinates();
                
                const Weights_T & r = C.EdgeLengths();
                
                const Int edge_count = C.EdgeCount();
                
                Vector_T b;
                
                // Handle wrap-around.
                for( Int i = 0; i < AmbDim; ++i )
                {
                    b[i] = r[edge_count-1] * ( p[edge_count][i]+p[0][i] );
                }
                
                for( Int k = 0; k < edge_count-1; ++ k )
                {
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        b[i] += r[i] * ( p[k][i] + p[k+1][i] );
                    }
                }
                
                const Real factor = static_cast<Real>(0.5)/r.Total();
                
                return b.Norm() * factor;
            }
            
            virtual Real MinValue( const Sampler_T & C ) const override
            {
                return static_cast<Real>(0);
            }
            
            virtual Real MaxValue( const Sampler_T & C ) const override
            {
                return Total(C.Omega());
            }
            
        public:
            
            virtual std::string Tag() const  override
            {
                return TO_STD_STRING(CLASS);
            }
        };
        
#undef CLASS
        
} // namespace CycleSampler
