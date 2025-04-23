#include "CoBarS.hpp"

// This demonstrate how to create an instance of `CoBarS::Sampler` and how to call its main method `RandomClosedPolygons` to generate random polygons along with their weights.


int main()
{
    // Dimensions of the ambient space has to be a compile-time constant.
    // Must be a compile-time constant
    constexpr int d = 3;
    
    const std::size_t edge_count = 64;
    
    // Create an instance of the CoBarS sampler.
    // This assumes that equilateral polygons shall be created.
    CoBarS::Sampler<d,double,std::size_t> S (edge_count);
    
    // CoBarS::Sampler has state and is thread-safe.
    // When you parallelize its use, then you should create an instance on each thread.
    // However, we see here an example that automatically parallelizes the call to `RandomClosedPolygons`.

    const std::size_t sample_count = 10000000;

    // Create containers for the polygons `p` and sampling weights `K`.
    std::vector<double> p ( sample_count * (edge_count + 1) * d );
    std::vector<double> K ( sample_count );
    
    // Note that `CoBarS::Sampler` is a biased sampler, so that you really have to use the sampling weights `K`; choose the valued of `quot_space_Q` depending on which probability distribution you would like to use.
        
    // Whether we want the sampling weights for the quotient space of polygons modulo SO(d) (`quot_space_Q = true`) or not (`quot_space_Q = false`).
    const bool quot_space_Q = true;
    
    // The worker routine is `CreateRandomClosedPolygons`.
    // It writes to raw pointers, so you are free to use whatever container you like,
    // as long as it stores its entries contiguously.
    
    
    // Generate random polygons using `thread_count` threads.
    
    const std::size_t thread_count = 8;
    
    
    std::cout << "Sampling " << sample_count <<" equilateral " << edge_count << "-gons, using " << thread_count << " threads." << std::endl;
    
    Tools::TimeInterval T_equilateral;
    
    T_equilateral.Tic();
    
    S.CreateRandomClosedPolygons(
        &p[0], &K[0], sample_count, quot_space_Q, thread_count
    );
    
    T_equilateral.Toc();
    
    std::cout << "Sampling done. Time elapsed = " << T_equilateral.Duration() << "." << std::endl;
    
    
    
    // Let's compute the sample average of the squared gyradius.
    {
        double weighted_sum = 0;
        double total_K      = 0;
        
        for( std::size_t k = 0; k < sample_count; ++k )
        {
            // Computing the squared gyradius.
            // CoBarS produces polygons with vanishing center of mass.
            // So the squared gyradius of a polygon (p_1,...,p_n) is given by
            // r^2 = \frac{1}{n} \sum_{i=1}^{n} |p_i|^2
            double r2 = 0;
            
            // Pointer to first entry of first vertex position of k-th polygon.
            const double * const p_ = &p[ (edge_count + 1) * d * k ];
            
            // We must add a contribution for the last vertex, because it is
            // just a copy of the first one.
            for( std::size_t i = 0; i < edge_count; ++i )
            {
                // Pointer to first entry of i-th vertex position.
                const double * const p_i = &p_[d * i];
                
                for( std::size_t j = 0; j < d; ++j )
                {
                    r2 += p_i[j] * p_i[j];
                }
            }
            
            r2 /= edge_count;
            
            weighted_sum += K[k] * r2;
            total_K      += K[k];
        }
        
        double weighted_mean = weighted_sum / total_K;
        
        std::cout << "Sample average of squared gyradius = " << weighted_mean << "." << std::endl;
        
        std::cout << "Theoretical expectation of squared gyradius = " <<  (edge_count+1.)/12. << "." << std::endl;
    }
    
    std::cout << std::endl;
    
    
    // Next we create an instance of the CoBarS sampler for nonequilateral polygons.
    // In addition to the `edge_count`, CoBarS::Sampler needs to know:
    
    // the edge lengths (stored here in the array `r`),...
    std::vector<double> r   ( edge_count );
    // ... and the weights for the Riemannian metric (stored here in the array `rho`).
    std::vector<double> rho ( edge_count );
    
    r[0]   = 2.;
    rho[0] = 1.;
    
    for( std::size_t i = 1; i < edge_count; ++i )
    {
        r[i]   = 1.;
        rho[i] = 1.;
    }
    
    CoBarS::Sampler<d,double,std::size_t> T ( &r[0], &rho[0], edge_count );
    
    std::cout << "Sampling " << sample_count <<" nonequilateral " << edge_count << "-gons, using " << thread_count << " threads." << std::endl;
    
    // The worker routine is `CreateRandomClosedPolygons`.
    // It writes to raw pointers, so you are free to use whatever container you like,
    // as long as it stores its entries contiguously.
    
    Tools::TimeInterval T_nonequilateral;
    
    T_nonequilateral.Tic();
    
    // This call is automatically parallelized over `thread_count` threads.
    T.CreateRandomClosedPolygons(
        &p[0], &K[0], sample_count, quot_space_Q, thread_count
    );
    
    T_nonequilateral.Toc();
    
    std::cout << "Sampling done. Time elapsed = " << T_nonequilateral.Duration() << "." << std::endl;
    
    // You can also sample a few preimplemented random variables on the polygon space without exporting the polygons themselves. This may safe a lot of memory traffic.
    
    // See Example_Sample_Binned/main.cpp for such an example.
    
    
    return 0;
}

