#pragma once
#include <cstdint>
#include <array>
#include <limits>
#include <type_traits>


// This is just a thin wrapper for std::mt19937_64 with the constructor that properly seeds its state from std::random_device;

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
    
    friend bool operator ==(const PCG64& lhs, const PCG64& rhs) noexcept
    {
        return (lhs.random_engine == rhs.random_engine);
    }
    
    friend bool operator !=(const PCG64& lhs, const PCG64& rhs) noexcept
    {
        return (lhs.random_engine != rhs.random_engine);
    }
    
    std::string ClassName()
    {
        return std::string("PCG64");
    }
};
