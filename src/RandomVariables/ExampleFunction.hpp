#pragma once

template<int AmbDim, typename Real = double, typename Int = long long>
class ExampleFunction : public RandomVariable<AmbDim,Real,Int>
{
public:
    
<<<<<<< HEAD
    // Imports several type aliases from the parent object.
    using Sampler_T         = typename BASE::Sampler_T;
    using SpherePoints_T    = typename BASE::SpherePoints_T;
    using SpacePoints_T     = typename BASE::SpacePoints_T;
    using Weights_T         = typename BASE::Weights_T;
    using Vector_T          = typename BASE::Vector_T;
=======
    using Sampler_T = Sampler<AmbDim,Real,Int>;
>>>>>>> 669f74e1da2608282dcd7df5c05e033802e4cfa6
    
    // Use this default constructor or write your own one.
    ExampleFunction() = default;
    
    // Use this default destructor or write your own one.
    virtual ~ExampleFunction() = default;
    
    // This inserts code for the Clone() routine.
    __ADD_CLONE_CODE__(ExampleFunction)

protected:
    
    // This is almost the only routine that you have to modify.
    // Calling the class on an object of class Sampler_T, it can access all its public data.
    // We provide hooks for the most interesting of them.
    virtual Real operator()( const Sampler_T & C ) const override
    {
        const Int n = C.EdgeCount();
        
        // The space coordinates as a (n+1) x AmbDim matrix.
        const Real * restrict const p = C.SpaceCoordinates();
        
        // The edge lengths as a vector of length n.
        const Real * restrict const edge_lengths = C.EdgeLengths();
        
        // The unit edge vectors of the conformally barycentered curve (array of size n x AmbDim).
        const Real * restrict const y = C.EdgeCoordinates();

        // The unit edge vectors of the open polygonal curve (array of size n x AmbDim).
        const Real * restrict const x = C.InitialEdgeCoordinates();

        // The shift vector (conformal barycenter of x w.r.t. edge_lengths) (array of size n AmbDim).
        const Real * restrict const w = C.ShiftVector();

        // DO SOMETHING MEANINGFUL HERE.
        
        return static_cast<Real>(0);
    }
    
    // Optionally, you can provide a lower bound for the range; this migh help with binning.
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
        return "ExampleFunction";
    }
};
