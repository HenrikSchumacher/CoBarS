#include "CoBarS.hpp"

// This demonstrate how to create an instance of `CoBarS::Sampler` and how to call its main method `RandomClosedPolygons` to generate random polygons along with their weights.


int main(int argc, const char * argv[])
{
    // Choose your base data types.
    // We recommend to use double and int, but you are free yo pick others.
    using Real = double;
    using Int  = int;

    // Dimensions of the ambient space has to be a compile-time constant.
    // Must be a compile-time constant
    constexpr Int d = 3;
    
    const     Int edge_count        = 64;
    const     Int sample_count      = 10000000;
    const     Int thread_count      = 8;
    
    // Create an instance of the CoBarS sampler.
    // This assumes that equilateral polygons shall be created.
    CoBarS::Sampler<d,Real,Int> S (edge_count);
    
    // CoBarS::Sampler has state and is thread-safe.
    // When you parallelize its use, then you should create an instance on each thread.
    // However, we see here an example that automatically parallelizes the call to `RandomClosedPolygons`.

    // Create containers for the data samples, more precisely for...
    
    // ... the edge vectors of open polygons, ...
    // ( to be interpreted as 3-tensor of size sample_count x edge_count x d  );
    std::vector<Real> x ( sample_count * edge_count * d );
    
    // ... the conformal barycenters, ...
    // ( to be interpreted as 2-tensor of size sample_count x d  );
    std::vector<Real> w ( sample_count * d );
    
    // ... the unit edge vectors of closed polygons, ...
    // ( to be interpreted as 3-tensor of size sample_count x edge_count x d  );
    std::vector<Real> y ( sample_count * edge_count * d );
    
    // ... the sample weights for the Pol space, ...
    std::vector<Real> K ( sample_count );
    
    // ... and the sample weights for the quotient space.
    std::vector<Real> K_quot ( sample_count );
    
    // Note that `CoBarS::Sampler` is a biased sampler, so that you really have to use the weights `K` or `K_quot`, depending on which probability distribution you would like to use.
    
    std::cout << "Sampling " << sample_count <<" equilateral " << edge_count << "-gons, using " << thread_count << " threads." << std::endl;
    
    // The worker routine is `RandomClosedPolygons`.
    // It writes to raw pointers, so you are free to use whatever container you like,
    // as long as it stores its entries contiguously.
    
    // This call is automatically parallelized over `thread_count` threads.
    S.RandomClosedPolygons(
        &x[0], &w[0], &y[0], &K[0], &K_quot[0], sample_count, thread_count
    );
    
    std::cout << "Done." << std::endl;

    
    
    
    
    
    
    // Create an instance of the CoBarS sampler for nonequilateral polygons.
    // In addition to the `edge_count`, CoBarS::Sampler needs to know...
    
    // ... the edge lengths (stored here in the array `r`),...
    std::vector<Real> r   ( edge_count );
    // ... the weights for the Riemannian metric (stored here in the array `rho`).
    std::vector<Real> rho ( edge_count );
    
    
    for( Int i = 0; i < edge_count; ++i )
    {
        r  [i] = 1 + (i/2);
        rho[i] = 1 + (i/2);
    }
    
    CoBarS::Sampler<d,Real,Int> T ( &r[0], &rho[0], edge_count );
    
    std::cout << "Sampling " << sample_count <<" nonequilateral " << edge_count << "-gons, using " << thread_count << " threads." << std::endl;
    
    // The worker routine is `RandomClosedPolygons`.
    // It writes to raw pointers, so you are free to use whatever container you like,
    // as long as it stores its entries contiguously.
    
    // This call is automatically parallelized over `thread_count` threads.
    T.RandomClosedPolygons(
        &x[0], &w[0], &y[0], &K[0], &K_quot[0], sample_count, thread_count
    );
    
    std::cout << "Done." << std::endl;
    
    
    
    // You can also sample a few preimplemented random variables on the polygon space without exporting the polygons themselves. This may safe a lot of memory traffic.
    
    // See Example_Sample_Binned/main.cpp for such an example.
    
    return 0;
}

