#pragma once

#include "../deps/pcg-cpp/include/pcg_random.hpp"


// This is just a thin wrapper for std::mt19937_64 with a constructor that properly seeds its state from std::random_device.

namespace CoBarS
{
    class PCG64
    {
    private:
        
        pcg64 random_engine;
        
    public:
        
        using result_type = pcg64::result_type;
        using UInt        = pcg64::result_type;
        
        PCG64() noexcept
        :   random_engine( pcg_extras::seed_seq_from<std::random_device>() )
        {}
        
        result_type operator()() noexcept
        {
            return random_engine();
        }
        
        static constexpr result_type min() noexcept
        {
            return pcg64::min();
        }
        
        static constexpr result_type max() noexcept
        {
            return pcg64::max();
        }
        
        std::string ClassName()
        {
            return std::string("PCG64");
        }
    };
    
} // namespace CoBarS
