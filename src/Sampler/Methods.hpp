public:

    using Base_T::Settings;
      
    Int EdgeCount() const override
    {
      return edge_count;
    }

    virtual const Weights_T & EdgeLengths() const override
    {
        return r;
    }
    
    virtual void ReadEdgeLengths( const Real * const r_in ) override
    {
        r.Read(r_in);
        
        total_r_inv = Inv( r.Total() );
    }
    
    
    virtual const Weights_T & Rho() const override
    {
        return rho;
    }
    
    virtual void ReadRho( const Real * const rho_in ) override
    {
        rho.Read(rho_in);
    }
    
    
    virtual void ComputeShiftVector() override
    {
        w.SetZero();
        
        if constexpr ( zerofyfirstQ )
        {
            for( Int k = 0; k < edge_count; ++k )
            {
                Vector_T x_k ( x, k );
                
                const Real r_k = r[k];
                
                for( Int i = 0; i < AmbDim; ++i )
                {
                    w[i] += x_k[i] * r_k;
                }
            }
        }
        else
        {
            // Overwrite by first summand.
            {
                Vector_T x_k ( x, 0 );
                
                const Real r_k = r[0];

                for( Int i = 0; i < AmbDim; ++i )
                {
                    w[i] = x_k[i] * r_k;
                }
            }

            // Add-in the others.
            for( Int k = 1; k < edge_count; ++k )
            {
                Vector_T x_k ( x, k );
                
                const Real r_k = r[k];

                for( Int i = 0; i < AmbDim; ++i )
                {
                    w[i] += x_k[i] * r_k;
                }
            }
        }

        
        // Normalize in that case that r does not sum up to 1.
        w *= total_r_inv;
    }
    
    virtual void ReadShiftVector( const Real * const w_in ) override
    {
        w.Read(w_in);
        
        // Use Euclidean barycenter as initial guess if the supplied initial guess does not make sense.
        if( Dot(w,w) > small_one )
        {
            ComputeShiftVector();
        }
    }
    
    virtual void ReadShiftVector( const Real * const w_in, const Int k ) override
    {
        ReadShiftVector( &w_in[ AmbDim * k ] );
    }
    
    virtual void WriteShiftVector( Real * w_out ) const override
    {
        w.Write( w_out );
    }
    
    virtual void WriteShiftVector( Real * w_out, const Int k ) const override
    {
        w.Write( &w_out[ AmbDim * k] );
    }
    
    virtual const Vector_T & ShiftVector() const override
    {
        return w;
    }
    
    
    virtual Real Residual() const override
    {
        return residual;
    }
    
    virtual Real ErrorEstimator() const override
    {
        return errorestimator;
    }
    
    virtual Int IterationCount() const override
    {
        return iter;
    }
    
    virtual Int MaxIterationCount() const override
    {
        return Settings().max_iter;
    }


private:

    static Real tanhc( const Real t )
    {
        // Computes tanh(t)/t in a stable way by using a Padé approximation around t = 0.
        constexpr Real a0 = Scalar::One<Real>;
        constexpr Real a1 = Frac<Real>(7,51);
        constexpr Real a2 = Frac<Real>(1,255);
        constexpr Real a3 = Frac<Real>(2,69615);
        constexpr Real a4 = Frac<Real>(1,34459425);
        
        constexpr Real b0 = Scalar::One<Real>;
        constexpr Real b1 = Frac<Real>(8,17);
        constexpr Real b2 = Frac<Real>(7,255);
        constexpr Real b3 = Frac<Real>(4,9945);
        constexpr Real b4 = Frac<Real>(1,765765);
        
        const Real t2 = t * t;
        
        const Real result = ( t2 <= one )
        ? (
            a0 + t2 * (a1 + t2 * (a2 + t2 * (a3 + t2 * a4)))
        )/(
            b0 + t2 * (b1 + t2 * (b2 + t2 * (b3 + t2 * b4)))
        )
        : ( t2 <= static_cast<Real>(7) ) ? std::tanh(t)/t : one/std::abs(t);
        
        return result;
    }


public:

    void NormalizeCompressedSamples(
        mut<Real> bins,
        const Int bin_count,
        mut<Real> moments_,
        const Int moment_count,
        const Int fun_count
    ) const
    {
        ptic(ClassName()+"::NormalizeCompressedSamples");
        for( Int i = 0; i < 3; ++i )
        {
            for( Int j = 0; j < fun_count; ++j )
            {
                // Normalize bins and moments.
                
                mut<Real> bins_i_j = &bins[ (i*fun_count+j) * bin_count ];
                
                mut<Real> moments_i_j = &moments_[ (i*fun_count+j) * moment_count ];
                
                // The field for zeroth moment is assumed to contain the total mass.
                Real factor = Inv( moments_i_j[0] );
                
                scale_buffer( factor, bins_i_j,    bin_count    );
                
                scale_buffer( factor, moments_i_j, moment_count );
            }
        }
        ptoc(ClassName()+"::NormalizeCompressedSamples");
    }
