#pragma once


namespace CoBarS
{
    
    // Just fill in the name of your new class in the next line; it will be automatically filled in below by the preprocessor.
    template<typename SamplerBase_T> class ExampleFunction;
    
    template<int AmbDim, typename Real, typename Int>
    class ExampleFunction<SamplerBase<AmbDim,Real,Int>>
    :   public RandomVariable<SamplerBase<AmbDim,Real,Int>>
    {
            
    public:
        
        using SamplerBase_T     = SamplerBase<AmbDim,Real,Int>;
        
    private:
        
        using Base_T            = RandomVariable<SamplerBase_T>;
        
    public:
        
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        // Use this default constructor or write your own one.
        ExampleFunction()
        {}
        
        // Use this default destructor or write your own one.
        virtual ~ExampleFunction()
        {}

    public:
        
        [[nodiscard]] std::shared_ptr<ExampleFunction> Clone () const
        {
            return std::shared_ptr<ExampleFunction>(CloneImplementation());
        }
                                                                                    
    private:
        
        [[nodiscard]] virtual ExampleFunction * CloneImplementation() const override
        {
            return new ExampleFunction(*this);
        }
        
    protected:
        
        // This is almost the only routine that you have to modify.
        // Calling the class on an object of class SamplerBase_T, it can access all its public data.
        // We provide hooks for the most interesting of them.
        virtual Real operator()( const SamplerBase_T & C ) const override
        {
            const Int edge_count = C.EdgeCount();
            
            // The edge lengths as a vector of length edge_count.
            const Weights_T & edge_lengths = C.EdgeLengths();
            
            // The shift vector (conformal barycenter of x w.r.t. edge_lengths).
            const Vector_T & w = C.ShiftVector();
            
            // DO SOMETHING MEANINGFUL HERE.
            
            return 0;
        }
        
        // Optionally, you can provide a lower bound for the range; this might help with binning.
        virtual Real MinValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 0;
        }
        
        // Optionally, you can provide an upper bound for the range; this migh help with binning.
        virtual Real MaxValue( const SamplerBase_T & C ) const override
        {
            (void)C;
            
            return 1;
        }
        
    public:
        
        // The name of the function that will appear in generated statistics. Don't modify this; it will be filled-in by the preprocessor.
        virtual std::string Tag() const  override
        {
            return std::string("ExampleFunction");
        }
    };
    
}  // namespace CoBarS
