#pragma once


namespace CycleSampler
{
    
    // Just fill in the name of your new class in the next line; it will be automatically filled in below by the preprocessor.
#define CLASS ExampleFunction
    
    template<int AmbDim, typename Real = double, typename Int = long long>
    class CLASS : public RandomVariable<AmbDim,Real,Int>
    {
    private:
        
        // It's not a bad thing to have a symbol for the base class...
        using Base_T            = RandomVariable<AmbDim,Real,Int>;
        
    public:
        
        // Imports several type aliases from the parent object.
        
        using Sampler_T         = typename Base_T::Sampler_T;
        using SpherePoints_T    = typename Base_T::SpherePoints_T;
        using SpacePoints_T     = typename Base_T::SpacePoints_T;
        using Weights_T         = typename Base_T::Weights_T;
        using Vector_T          = typename Base_T::Vector_T;
        
        // Use this default constructor or write your own one.
        CLASS()
        {}
        
        // Use this default destructor or write your own one.
        virtual ~CLASS()
        {}
        
        // This inserts code for the Clone() routine.
        __ADD_CLONE_CODE__(CLASS)
        
    protected:
        
        // This is almost the only routine that you have to modify.
        // Calling the class on an object of class Sampler_T, it can access all its public data.
        // We provide hooks for the most interesting of them.
        virtual Real operator()( const Sampler_T & C ) const override
        {
            const Int edge_count = C.EdgeCount();
            
            // The space coordinates as a edge_count x AmbDim matrix.
            const SpacePoints_T & p = C.SpaceCoordinates();
            
            // The edge lengths as a vector of length edge_count.
            const Weights_T & edge_lengths = C.EdgeLengths();
            
            // The unit edge vectors of the conformally barycentered curve.
            const SpherePoints_T & y = C.EdgeCoordinates();
            
            // The unit edge vectors of the open polygonal curve.
            const SpherePoints_T & x = C.InitialEdgeCoordinates();
            
            // The shift vector (conformal barycenter of x w.r.t. edge_lengths).
            const Vector_T & w = C.ShiftVector();
            
            // DO SOMETHING MEANINGFUL HERE.
            
            return static_cast<Real>(0);
        }
        
        // Optionally, you can provide a lower bound for the range; this might help with binning.
        virtual Real MinValue( const Sampler_T & C ) const override
        {
            return static_cast<Real>(0);
        }
        
        // Optionally, you can provide an upper bound for the range; this migh help with binning.
        virtual Real MaxValue( const Sampler_T & C ) const override
        {
            return static_cast<Real>(1);
        }
        
    public:
        
        // The name of the function that will appear in generated statistics. Don't modify this; it will be filled-in by the preprocessor.
        virtual std::string Tag() const  override
        {
            return TO_STD_STRING(CLASS);
        }
    };
    
#undef CLASS
    
}
