public:

    using Base_T::Settings;
      
    Int EdgeCount() const override
    {
      return edge_count_;
    }

    virtual const Weights_T & EdgeLengths() const override
    {
        return r_;
    }
    
    virtual void ReadEdgeLengths( const Real * const r ) override
    {
        r_.Read(r);
        
        total_r_inv = Inv( r_.Total() );
    }
    
    
    virtual const Weights_T & Rho() const override
    {
        return rho_;
    }
    
    virtual void ReadRho( const Real * const rho ) override
    {
        rho_.Read(rho);
    }
    
private:
    
    virtual void ComputeInitialShiftVector() override
    {
        if constexpr ( zerofyfirstQ )
        {
            w_.SetZero();
            
            for( Int i = 0; i < edge_count_; ++i )
            {
                Vector_T x_i ( x_, i );
                
                const Real r_i = r_[i];
                
                for( Int j = 0; j < AmbDim; ++j )
                {
                    w_[j] += x_i[j] * r_i;
                }
            }
        }
        else
        {
            // Overwrite by first summand.
            {
                Vector_T x_i ( x_, 0 );
                
                const Real r_i = r_[0];

                for( Int j = 0; j < AmbDim; ++j )
                {
                    w_[j] = x_i[j] * r_i;
                }
            }

            // Add-in the others.
            for( Int i = 1; i < edge_count_; ++i )
            {
                Vector_T x_i ( x_, i );
                
                const Real r_i = r_[i];

                for( Int j = 0; j < AmbDim; ++j )
                {
                    w_[j] += x_i[j] * r_i;
                }
            }
        }

        
        // Normalize in that case that r does not sum up to 1.
        w_ *= total_r_inv;
    }

public:
    
    virtual void ReadShiftVector( const Real * restrict const w, const Int offset = 0 ) override
    {
        w_.Read( &w[ AmbDim * offset ] );
        
        // Use Euclidean barycenter as initial guess if the supplied initial guess does not make sense.
        if( Dot(w_,w_) > small_one )
        {
            ComputeInitialShiftVector();
        }
    }
    
    virtual void WriteShiftVector( Real * restrict w, const Int offset = 0 ) const override
    {
        w_.Write( &w[ AmbDim * offset] );
    }
    
    virtual const Vector_T & ShiftVector() const override
    {
        return w_;
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
        constexpr Real a0 = one;
        constexpr Real a1 = Frac<Real>(7,51);
        constexpr Real a2 = Frac<Real>(1,255);
        constexpr Real a3 = Frac<Real>(2,69615);
        constexpr Real a4 = Frac<Real>(1,34459425);
        
        constexpr Real b0 = one;
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
        mptr<Real> bins, const Int bin_count,
        mptr<Real> moms, const Int mom_count,
        const Int fun_count
    ) const
    {
        ptic(ClassName()+"::NormalizeCompressedSamples");
        for( Int i = 0; i < 3; ++i )
        {
            for( Int j = 0; j < fun_count; ++j )
            {
                // Normalize bins and moments.
                
                mptr<Real> bins_i_j = &bins[ (i*fun_count+j) * bin_count ];
                
                mptr<Real> moments_i_j = &moms[ (i*fun_count+j) * mom_count ];
                
                // The field for zeroth moment is assumed to contain the total mass.
                Real factor = Inv( moments_i_j[0] );
                
                scale_buffer( factor, bins_i_j,    bin_count    );
                
                scale_buffer( factor, moments_i_j, mom_count );
            }
        }
        ptoc(ClassName()+"::NormalizeCompressedSamples");
    }
