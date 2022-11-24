#pragma once

namespace CycleSampler {

    template<int AmbDim, typename Real = double, typename Int = long long>
    class ShiftMap
    {
        ASSERT_FLOAT(Real);
        ASSERT_INT(Int);
        
    public:
        
        ShiftMap() = default;
        
        virtual ~ShiftMap(){}
        
        static constexpr Real zero              = 0;
        static constexpr Real half              = 0.5;
        static constexpr Real one               = 1;
        static constexpr Real two               = 2;
        static constexpr Real three             = 3;
        static constexpr Real four              = 4;
        static constexpr Real eps               = std::numeric_limits<Real>::min();
        static constexpr Real infty             = std::numeric_limits<Real>::max();
        static constexpr Real small_one         = 1 - 16 * eps;
        static constexpr Real big_one           = 1 + 16 * eps;
        static constexpr Real g_factor          = 4;
        static constexpr Real g_factor_inv      = one/g_factor;
        static constexpr Real norm_threshold    = 0.99 * 0.99 + 16 * eps;
        static constexpr Real two_pi            = static_cast<Real>(2 * M_PI);
        
    protected:
        
        SmallSquareMatrix<AmbDim,Real,Int> cbar_mat;
        SmallSquareMatrix<AmbDim,Real,Int> gamma_mat;
        
    public:
        
        using Vector_T          = SmallVector<AmbDim,Real,Int>;
        using SquareMatrix_T    = SmallSquareMatrix<AmbDim,Real,Int>;
        using SymmetricMatrix_T = SmallSymmetricMatrix<AmbDim,Real,Int>;
        
        void Shift(
            const Real * restrict const x_in,
            const Real * restrict const s,
                  Real * restrict const y_out,
            const Int                   n,
            const Real                  t = one
        ) const
        {
            // Shifts all entries of x along y and writes the results to y.
            // Mind that x and y are stored in SoA fashion, i.e., as matrix of size AmbDim x point_count.
            
            Real s_[AmbDim];
            Real z_[AmbDim];
                        
            Real ss = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                s_[i] = t * s[i];
                
                ss += s_[i] * s_[i];
            }
            
            const Real one_minus_ss = big_one - ss;
            const Real one_plus_ss  = big_one + ss;

            if( ss <= norm_threshold )
            {
//                #pragma clang loop vectorize(assume_safety)
                for( Int k = 0; k < n; ++k )
                {
                    Real sz2 = 0;
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        z_[i] = x_in[AmbDim*k+i];
                        
                        sz2 += s_[i] * z_[i];
                    }
                    
                    sz2 *= two;
                    
                    const Real denom = one / ( one_plus_ss - sz2 );
                    
                    const Real sz2_minus_2 = sz2 - two;

                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        y_out[AmbDim*k+i] = (one_minus_ss * z_[i] + sz2_minus_2 * s_[i]) * denom;
                    }
                }
            }
            else
            {
                // If w lies close to the boundary of the ball, then normalizing the output is a good idea.

//                #pragma clang loop vectorize(assume_safety)
                for( Int k = 0; k < n; ++k )
                {
                    Real sz2 = 0;
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        z_[i] = x_in[AmbDim*k+i];
                        
                        sz2 += s_[i] * z_[i];
                    }
                    
                    sz2 *= two;
                    
                    const Real denom = one / ( one_plus_ss - sz2 );
                    
                    const Real sz2_minus_2 = sz2 - two;

                    Real r2 = 0;
                    
                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        z_[i] = (one_minus_ss * z_[i] + sz2_minus_2 * s_[i]) * denom;
                        
                        r2 += z_[i] * z_[i];
                    }
                    
                    const Real z_inv_norm = one/sqrt(r2);

                    for( Int i = 0; i < AmbDim; ++i )
                    {
                        y_out[AmbDim*k+i] = z_[i] * z_inv_norm;
                    }
                }
            }
        }
        
        void DShiftD1(
            const Real * restrict const x_in,
            const Real * restrict const s,
                  Real * restrict const a_out,
            const Int                   n,
            const Real                  t = one
        ) const
        {
            // Shifts all entries of x along y and writes the results to y.
            // Mind that x and y are stored in SoA fashion, i.e., as matrix of size AmbDim x point_count.
            
            Real s_[AmbDim];
            Real z_[AmbDim];
            
            Real ss = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                s_[i] = t * s[i];
                
                ss += s_[i] * s_[i];
            }
            
            const Real one_minus_ss = big_one - ss;
            const Real one_plus_ss  = big_one + ss;
           
            for( Int k = 0; k < n; ++k )
            {
                Real sz2 = 0;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    z_[i] = x_in[AmbDim*k+i];
                    
                    sz2 += s_[i] * z_[i];
                }
                
                sz2 *= two;
                
                const Real denom = one / ( one_plus_ss - sz2 );
                
                const Real sz2_minus_2 = sz2 - two;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    for( Int j = 0; j < AmbDim; ++j )
                    {
                        a_out[ AmbDim*AmbDim * k + AmbDim * i + j ] = t * two * denom * (
                            s_[i] * z_[j] - z_[i] * s_[j]
                            +
                            denom * ( one_minus_ss * z_[i] + sz2_minus_2 * s_[i]) * (z_[j] - s_[j])
                        );
                    }
                }
                
                const Real diag = t * sz2_minus_2 * denom;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    a_out[ AmbDim*AmbDim * k + AmbDim * i + i ] += diag;
                }
            }
        }
        
        
        void DShiftD2(
            const Real * restrict const x_in,
            const Real * restrict const s,
                  Real * restrict const b_out,
            const Int                   n,
            const Real                  t = one
        ) const
        {
            // Shifts all entries of x along y and writes the results to y.
            // Mind that x and y are stored in SoA fashion, i.e., as matrix of size AmbDim x point_count.
            
            Real s_[AmbDim];
            Real z_[AmbDim];
            
            Real ss = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                s_[i] = t * s[i];
                
                ss += s_[i] * s_[i];
            }
            
            const Real one_minus_ss = big_one - ss;
            const Real one_plus_ss  = big_one + ss;
           
            for( Int k = 0; k < n; ++k )
            {
                Real sz2 = 0;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    z_[i] = x_in[AmbDim*k+i];
                    
                    sz2 += s_[i] * z_[i];
                }
                
                sz2 *= two;
                
                const Real denom = one / ( one_plus_ss - sz2 );
                
                const Real sz2_minus_2 = sz2 - two;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    for( Int j = 0; j < AmbDim; ++j )
                    {
                        b_out[ AmbDim*AmbDim * k + AmbDim * i + j ] = t * two * denom * (
                            - s_[i] * (z_[j]-s_[j])
                            +
                            denom * ( one_minus_ss * z_[i] + sz2_minus_2 * s_[i]) * (s_[j] - ss * z_[j])
                        );
                    }
                }
                
                const Real diag  = t * one_minus_ss * denom;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    b_out[ AmbDim*AmbDim * k + AmbDim * i + i ] += diag;
                }
            }
        }
        
        void Compute_c(
            const Real * restrict const y_in,
            const Real * restrict const w_in,
            const Real * restrict const x_in,
                  Real * restrict const c_out,
            const Int                   n
        ) const
        {
            // Shifts all entries of x along y and writes the results to y.
            // Mind that x and y are stored in SoA fashion, i.e., as matrix of size AmbDim x point_count.
            
            Real w [AmbDim];
            
            Real ww = 0;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                w[i]  = w_in[i];
                
                ww += w[i] * w[i];
            }
            
            const Real one_minus_ww = big_one - ww;
            const Real one_plus_ww  = big_one + ww;
           
            for( Int k = 0; k < n; ++k )
            {
                Real y    [AmbDim];
                Real x    [AmbDim];
                
                Real a    [AmbDim][AmbDim];
                Real b_inv[AmbDim][AmbDim];
                
                Real wy = 0;
                Real wx = 0;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    y[i] = y_in[AmbDim*k+i];
                    x[i] = x_in[AmbDim*k+i];
                    
                    wy += w[i] * y[i];
                    wx += w[i] * x[i];
                }
                
                const Real a_denom     = one / ( one_plus_ww + two * wy );
                const Real b_inv_denom = one / ( one_plus_ww - two * wx );
                
                
                const Real a_diag = two * (one + wy) * a_denom;
                const Real b_inv_diag = one_minus_ww * b_inv_denom;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    for( Int j = 0; j < AmbDim; ++j )
                    {
                        a[i][j] = 2 * a_denom * (
                            w[i] * y[j] - y[i] * w[j] - x[i] * ( y[j] + w[j] )
                        );
                        
                        b_inv[i][j] = two * b_inv_denom * (
                            y[i] * (w[j] - ww * x[j])- w[i] * (x[j]-w[j])
                        );
                    }
                    
                    a    [i][i] += a_diag;
                    b_inv[i][i] += b_inv_diag;
                }

                Real c    [AmbDim][AmbDim] = {};
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    for( Int l = 0; l < AmbDim; ++l )
                    {
                        for( Int j = 0; j < AmbDim; ++j )
                        {
                            c[i][j] += b_inv[i][l] * a[l][j];
                        }
                    }
                }
                
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    for( Int j = 0; j < AmbDim; ++j )
                    {
                        c_out[ AmbDim*AmbDim * k + AmbDim * i + j ] = c[i][j];
                    }
                }
            }
        }
        
        Real EdgeSpaceSamplingWeight(
            const Real * restrict const x_in,
            const Real * restrict const w_in,
            const Real * restrict const y_in,
            const Real * restrict const omega_in,
            const Real * restrict const rho_in,
            const Int                   n
        )
        {
            // Shifts all entries of x along y and writes the results to y.
            // Mind that x and y are stored in SoA fashion, i.e., as matrix of size AmbDim x point_count.
            
            Real cbar  [AmbDim][AmbDim] = {{}};
            Real gamma [AmbDim][AmbDim] = {{}};
            
            Real prod = one;
            
            Real w [AmbDim];
            
            Real ww = zero;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                w[i] = w_in[i];
                
                ww += w[i] * w[i];
            }
            
            //            const Real one_minus_ww    = big_one - ww;
            const Real one_plus_ww     = big_one + ww;
            //            const Real one_plus_ww_inv = one / one_plus_ww;
            
            for( Int k = 0; k < n; ++k )
            {
                Real y [AmbDim];
                
                Real wy = zero;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    y[i] = y_in[AmbDim*k+i];
                    
                    wy += w[i] * y[i];
                }
                
                const Real factor = one_plus_ww + two * wy;
                
                // Multiplying by one_plus_ww_inv so that prod does not grow so quickly.
                //                prod *= factor * one_plus_ww_inv;
                prod *= factor;
                
                const Real omega_k = omega_in[k];
                const Real omega_over_rho_k = omega_k/rho_in[k];
                const Real omega_over_rho_k_squared = omega_over_rho_k * omega_over_rho_k;
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    for( Int j = 0; j < AmbDim; ++j )
                    {
                        const Real scratch = (static_cast<Real>(i==j) - y[i]*y[j]);
                        
                        gamma[i][j] += omega_over_rho_k_squared * scratch;
                        
                        cbar[i][j]  += omega_k * scratch;
                    }
                }
            }
            
            gamma_mat.Read(&gamma[0][0]);
            cbar_mat.Read(&cbar[0][0]);
            
            return MyMath::pow( prod, static_cast<Int>(AmbDim - 1) ) * std::sqrt(gamma_mat.Det()) / cbar_mat.Det();
        }

    public:

        static constexpr Int AmbientDimension()
        {
            return AmbDim;
        }

        static std::string ClassName()
        {
            return "ShiftMap<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
        }
    };
    
} // namespace CycleSampler
  
//#pragma once
//
//namespace CycleSampler
//{
//
//    template<int AmbDim, typename Real = double, typename Int = long long>
//    class ShiftMap
//    {
//        ASSERT_FLOAT(Real);
//        ASSERT_INT(Int);
//
//    public:
//
//        ShiftMap() = default;
//
//        virtual ~ShiftMap() = default;
//
//        static constexpr Real zero              = 0;
//        static constexpr Real half              = 0.5;
//        static constexpr Real one               = 1;
//        static constexpr Real two               = 2;
//        static constexpr Real three             = 3;
//        static constexpr Real four              = 4;
//        static constexpr Real eps               = std::numeric_limits<Real>::min();
//        static constexpr Real infty             = std::numeric_limits<Real>::max();
//        static constexpr Real small_one         = 1 - 16 * eps;
//        static constexpr Real big_one           = 1 + 16 * eps;
//        static constexpr Real g_factor          = 4;
//        static constexpr Real g_factor_inv      = one/g_factor;
//        static constexpr Real norm_threshold    = 0.99 * 0.99 + 16 * eps;
//        static constexpr Real two_pi            = static_cast<Real>(2 * M_PI);
//
//    public:
//
//        void Shift(
//            const Real * restrict const x_in,
//            const Real * restrict const s,
//                  Real * restrict const y_out,
//            const Int                   n,
//            const Real                  t = one
//        ) const
//        {
//            // Shifts all entries of x along s and writes the results to y.
//
//            Real s_[AmbDim];
//            Real z_[AmbDim];
//
//            Real ss = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                s_[i] = t * s[i];
//
//                ss += s_[i] * s_[i];
//            }
//
//            const Real one_minus_ss = big_one - ss;
//            const Real one_plus_ss  = big_one + ss;
//
//            if( ss <= norm_threshold )
//            {
//                for( Int k = 0; k < n; ++k )
//                {
//                    Real sz2 = 0;
//
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        z_[i] = x_in[AmbDim*k+i];
//
//                        sz2 += s_[i] * z_[i];
//                    }
//
//                    sz2 *= two;
//
//                    const Real denom = one / ( one_plus_ss - sz2 );
//
//                    const Real sz2_minus_2 = sz2 - two;
//
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        y_out[AmbDim*k+i] = (one_minus_ss * z_[i] + sz2_minus_2 * s_[i]) * denom;
//                    }
//                }
//            }
//            else
//            {
//                // If w lies close to the boundary of the ball, then normalizing the output is a good idea.
//
//                for( Int k = 0; k < n; ++k )
//                {
//                    Real sz2 = 0;
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        z_[i] = x_in[AmbDim*k+i];
//
//                        sz2 += s_[i] * z_[i];
//                    }
//
//                    sz2 *= two;
//
//                    const Real denom = one / ( one_plus_ss - sz2 );
//
//                    const Real sz2_minus_2 = sz2 - two;
//
//                    Real r2 = 0;
//
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        z_[i] = (one_minus_ss * z_[i] + sz2_minus_2 * s_[i]) * denom;
//
//                        r2 += z_[i] * z_[i];
//                    }
//
//                    const Real z_inv_norm = one/sqrt(r2);
//
//                    for( Int i = 0; i < AmbDim; ++i )
//                    {
//                        y_out[AmbDim*k+i] = z_[i] * z_inv_norm;
//                    }
//                }
//            }
//        }
//
//        void DShiftD1(
//            const Real * restrict const x_in,
//            const Real * restrict const s,
//                  Real * restrict const a_out,
//            const Int                   n,
//            const Real                  t = one
//        ) const
//        {
//            Real s_[AmbDim];
//            Real z_[AmbDim];
//
//            Real ss = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                s_[i] = t * s[i];
//
//                ss += s_[i] * s_[i];
//            }
//
//            const Real one_minus_ss = big_one - ss;
//            const Real one_plus_ss  = big_one + ss;
//
//            for( Int k = 0; k < n; ++k )
//            {
//                Real sz2 = 0;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    z_[i] = x_in[AmbDim*k+i];
//
//                    sz2 += s_[i] * z_[i];
//                }
//
//                sz2 *= two;
//
//                const Real denom = one / ( one_plus_ss - sz2 );
//
//                const Real sz2_minus_2 = sz2 - two;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    for( Int j = 0; j < AmbDim; ++j )
//                    {
//                        a_out[ AmbDim*AmbDim * k + AmbDim * i + j ] = t * two * denom * (
//                            s_[i] * z_[j] - z_[i] * s_[j]
//                            +
//                            denom * ( one_minus_ss * z_[i] + sz2_minus_2 * s_[i]) * (z_[j] - s_[j])
//                        );
//                    }
//                }
//
//                const Real diag = t * sz2_minus_2 * denom;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    a_out[ AmbDim*AmbDim * k + AmbDim * i + i ] += diag;
//                }
//            }
//        }
//
//
//        void DShiftD2(
//            const Real * restrict const x_in,
//            const Real * restrict const s,
//                  Real * restrict const b_out,
//            const Int                   n,
//            const Real                  t = one
//        ) const
//        {
//            Real s_[AmbDim];
//            Real z_[AmbDim];
//
//            Real ss = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                s_[i] = t * s[i];
//
//                ss += s_[i] * s_[i];
//            }
//
//            const Real one_minus_ss = big_one - ss;
//            const Real one_plus_ss  = big_one + ss;
//
//            for( Int k = 0; k < n; ++k )
//            {
//                Real sz2 = 0;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    z_[i] = x_in[AmbDim*k+i];
//
//                    sz2 += s_[i] * z_[i];
//                }
//
//                sz2 *= two;
//
//                const Real denom = one / ( one_plus_ss - sz2 );
//
//                const Real sz2_minus_2 = sz2 - two;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    for( Int j = 0; j < AmbDim; ++j )
//                    {
//                        b_out[ AmbDim*AmbDim * k + AmbDim * i + j ] = t * two * denom * (
//                            - s_[i] * (z_[j]-s_[j])
//                            +
//                            denom * ( one_minus_ss * z_[i] + sz2_minus_2 * s_[i]) * (s_[j] - ss * z_[j])
//                        );
//                    }
//                }
//
//                const Real diag  = t * one_minus_ss * denom;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    b_out[ AmbDim*AmbDim * k + AmbDim * i + i ] += diag;
//                }
//            }
//        }
//
//        void Compute_c(
//            const Real * restrict const y_in,
//            const Real * restrict const w_in,
//            const Real * restrict const x_in,
//                  Real * restrict const c_out,
//            const Int                   n
//        ) const
//        {
//            Real w [AmbDim];
//
//            Real ww = 0;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                w[i] = w_in[i];
//
//                ww += w[i] * w[i];
//            }
//
//            const Real one_minus_ww = big_one - ww;
//            const Real one_plus_ww  = big_one + ww;
//
//            for( Int k = 0; k < n; ++k )
//            {
//                Real y     [AmbDim];
//                Real x     [AmbDim];
//
//                Real a     [AmbDim][AmbDim];
//                Real b_inv [AmbDim][AmbDim];
//
//                Real wy = 0;
//                Real wx = 0;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    y[i] = y_in[AmbDim*k+i];
//                    x[i] = x_in[AmbDim*k+i];
//
//                    wy += w[i] * y[i];
//                    wx += w[i] * x[i];
//                }
//
//                const Real a_denom     = one / ( one_plus_ww + two * wy );
//                const Real b_inv_denom = one / ( one_plus_ww - two * wx );
//
//
//                const Real a_diag = two * (one + wy) * a_denom;
//                const Real b_inv_diag = one_minus_ww * b_inv_denom;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    for( Int j = 0; j < AmbDim; ++j )
//                    {
//                        a[i][j] = 2 * a_denom * (
//                            w[i] * y[j] - y[i] * w[j] - x[i] * ( y[j] + w[j] )
//                        );
//
//                        b_inv[i][j] = two * b_inv_denom * (
//                            y[i] * (w[j] - ww * x[j])- w[i] * (x[j]-w[j])
//                        );
//                    }
//
//                    a    [i][i] += a_diag;
//                    b_inv[i][i] += b_inv_diag;
//                }
//
//                Real c    [AmbDim][AmbDim] = {};
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    for( Int l = 0; l < AmbDim; ++l )
//                    {
//                        for( Int j = 0; j < AmbDim; ++j )
//                        {
//                            c[i][j] += b_inv[i][l] * a[l][j];
//                        }
//                    }
//                }
//
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    for( Int j = 0; j < AmbDim; ++j )
//                    {
//                        c_out[ AmbDim*AmbDim * k + AmbDim * i + j ] = c[i][j];
//                    }
//                }
//            }
//        }
//
//        Real EdgeSpaceSamplingWeight(
//            const Real * restrict const x_in,
//            const Real * restrict const w_in,
//            const Real * restrict const y_in,
//            const Real * restrict const omega_in,
//            const Real * restrict const rho_in,
//            const Int                   n
//        ) const
//        {
//            Real cbar  [AmbDim][AmbDim] = {};
//            Real gamma [AmbDim][AmbDim] = {};
//
//            Real prod = one;
//
//            Real w [AmbDim];
//
//            Real ww = zero;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                w[i]  = w_in[i];
//
//                ww += w[i] * w[i];
//            }
//
////            const Real one_minus_ww    = big_one - ww;
//            const Real one_plus_ww     = big_one + ww;
////            const Real one_plus_ww_inv = one / one_plus_ww;
//
//            for( Int k = 0; k < n; ++k )
//            {
//                Real y [AmbDim];
//
//                Real wy = zero;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    y[i] = y_in[AmbDim*k+i];
//
//                    wy += w[i] * y[i];
//                }
//
//                const Real factor = one_plus_ww + two * wy;
//
//                // Multiplying by one_plus_ww_inv so that prod does not grow so quickly.
////                prod *= factor * one_plus_ww_inv;
//                prod *= factor;
//
//                const Real omega_k = omega_in[k];
//                const Real omega_over_rho_k = omega_k/rho_in[k];
//                const Real omega_over_rho_k_squared = omega_over_rho_k * omega_over_rho_k;
//
//                for( Int i = 0; i < AmbDim; ++i )
//                {
//                    for( Int j = 0; j < AmbDim; ++j )
//                    {
//                        const Real scratch = (static_cast<Real>(i==j) - y[i]*y[j]);
//
//                        gamma[i][j] += omega_over_rho_k_squared * scratch;
//
//                        cbar[i][j]  += omega_k * scratch;
//                    }
//                }
//            }
//
////             We can simply absorb the factor std::pow(2/(one_minus_ww),d) into the function chi.
////            {
////                const Real scratch = static_cast<Real>(2)/(one_minus_ww);
////
////                for( Int i = 0; i < AmbDim; ++i )
////                {
////                    for( Int j = 0; j < AmbDim; ++j )
////                    {
////                        cbar(i,j) *= scratch;
////                    }
////                }
////            }
//
//            SmallSquareMatrix<AmbDim,Real,Int> cbar_mat;
//            SmallSquareMatrix<AmbDim,Real,Int> gamma_mat;
//
//            for( Int i = 0; i < AmbDim; ++i )
//            {
//                for( Int j = 0; j < AmbDim; ++j )
//                {
//                    cbar_mat (i,j) = cbar [i][j];
//                    gamma_mat(i,j) = gamma[i][j];
//                }
//            }
//
//            return MyMath::pow( prod, static_cast<Int>(AmbDim - 1) ) * std::sqrt(gamma_mat.Det()) / cbar_mat.Det();
//        }
//
//    public:
//
//        static constexpr Int AmbientDimension()
//        {
//            return AmbDim;
//        }
//
//        static std::string ClassName()
//        {
//            return "ShiftMap<"+ToString(AmbDim)+","+TypeName<Real>::Get()+","+TypeName<Int>::Get()+">";
//        }
//    };
//
//} // namespace CycleSampler
