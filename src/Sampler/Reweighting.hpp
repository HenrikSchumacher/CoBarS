public:

    virtual Real EdgeSpaceSamplingWeight() const override
    {
        return edge_space_sampling_weight;
    }

    virtual Real EdgeQuotientSpaceSamplingCorrection() const override
    {
        return edge_quotient_space_sampling_correction;
    }

    virtual Real EdgeQuotientSpaceSamplingWeight() const override
    {
        return EdgeSpaceSamplingWeight() * EdgeQuotientSpaceSamplingCorrection();
    }

protected:

    void ComputeEdgeSpaceSamplingWeight()
    {
        // Shifts all entries of x along y and writes the results to y.
        // Mind that x and y are stored in SoA fashion, i.e., as matrix of size AmbDim x point_count.
        
        SquareMatrix_T cbar  (zero);
        SquareMatrix_T gamma (zero);
        
        Real prod = one;
        
        const Real ww = Dot(w,w);
        
        const Real one_plus_ww = big_one + ww;
        
        for( Int k = 0; k < edge_count; ++k )
        {
            Vector_T y_k ( y, k );
            
            const Real wy = Dot(w,y_k);
            
            const Real factor = one_plus_ww + two * wy;
            
            // Multiplying by one_plus_ww_inv so that prod does not grow so quickly.
            prod *= factor;
            
            const Real r_k = r[k];
            const Real r_over_rho_k = r_k / rho[k];
            const Real r_over_rho_k_squared = r_over_rho_k * r_over_rho_k;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                for( Int j = 0; j < AmbDim; ++j )
                {
                    const Real scratch = static_cast<Real>(i==j) - y_k[i] * y_k[j];
                    
                    gamma[i][j] += r_over_rho_k_squared * scratch;
                    
                    cbar [i][j] += r_k * scratch;
                }
            }
        }
        
        // We can simply absorb the factor std::pow(2/(one_minus_ww),d) into the function chi.
        //  cbar *= static_cast<Real>(2)/(one_minus_ww);
        
        edge_space_sampling_weight = Power(prod, static_cast<Int>(AmbDim-1)) * sqrt(gamma.Det()) / cbar.Det();
    }

    void ComputeEdgeQuotientSpaceSamplingCorrection()
    {
        if constexpr ( AmbDim == 2)
        {
            edge_quotient_space_sampling_correction = one;
            return;
        }
        
        Tiny::SelfAdjointMatrix<AmbDim, Real, Int> Sigma;
        
        // We fill only the upper triangle of Sigma, because that's the only thing that the function Eigenvalues needs.
      
        
        if constexpr ( zerofyfirstQ )
        {
            Sigma.SetZero();
            
            for( Int k = 0; k < edge_count; ++k )
            {
                Vector_T y_k ( y, k );
                
                const Real rho_squared = rho[k] * rho[k];
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = rho_squared * y_k[i];
                    
                    for( Int j = i; j < AmbDim; ++j )
                    {
                        Sigma[i][j] += factor * y_k[j];
                    }
                }
            }
        }
        else
        {
            // Overwrite for k = 0.
            {
                constexpr Int k = 0;
                
                Vector_T y_k ( y, k );
                
                const Real rho_squared = rho[0] * rho[0];
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = rho_squared * y_k[i];

                    for( Int j = i; j < AmbDim; ++j )
                    {
                        Sigma[i][j] = factor * y_k[j];
                    }
                }
            }

            // Now we add-in the other entries.
            for( Int k = 1; k < edge_count; ++k )
            {
                Vector_T y_k ( y, k );
                
                const Real rho_squared = rho[k] * rho[k];
                for( Int i = 0; i < AmbDim; ++i )
                {
                    const Real factor = rho_squared * y_k[i];

                    for( Int j = i; j < AmbDim; ++j )
                    {
                        Sigma[i][j] += factor * y_k[j];
                    }
                }
            }
        }
        
        
        
        if constexpr ( AmbDim == 3 )
        {
            // Exploiting that
            //      (lambda[0] + lambda[1]) * (lambda[0] + lambda[2]) * (lambda[1] + lambda[2])
            //      =
            //      ( tr(Sigma*Sigma) - tr(Sigma)*tr(Sigma) ) *  tr(Sigma)/2 - det(Sigma)
            //  Thus, it can be expressed by as third-order polynomial in the entries of the matrix.
            
            const Real S_00 = Sigma[0][0] * Sigma[0][0];
            const Real S_11 = Sigma[1][1] * Sigma[1][1];
            const Real S_22 = Sigma[2][2] * Sigma[2][2];
            
            const Real S_10 = Sigma[0][1] * Sigma[0][1];
            const Real S_20 = Sigma[0][2] * Sigma[0][2];
            const Real S_21 = Sigma[1][2] * Sigma[1][2];
            
            const Real det = std::abs(
                  Sigma[0][0] * ( S_11 + S_22 - S_10 - S_20 )
                + Sigma[1][1] * ( S_00 + S_22 - S_10 - S_21 )
                + Sigma[2][2] * ( S_00 + S_11 - S_20 - S_21 )
                + two * (Sigma[0][0]*Sigma[1][1]*Sigma[2][2] - Sigma[0][1]*Sigma[0][2]*Sigma[1][2])
            );
            edge_quotient_space_sampling_correction = one / std::sqrt(det);
        }
        else
        {
            Tiny::Vector<AmbDim,Real,Int> lambda;
            
            // Compute eigenvalues by QR algorithm.
            
            Sigma.Eigenvalues(lambda, Sqrt(Scalar::eps<Real>), 16);
            
            Real det = one;
            
            for( Int i = 0; i < AmbDim; ++i )
            {
                for( Int j = i+1; j < AmbDim; ++j )
                {
                    det *= (lambda(i)+lambda(j));
                }
            }
            
            edge_quotient_space_sampling_correction = one / std::sqrt(det);
        }
    }
