#pragma once
#include <cstdint>
#include <array>
#include <limits>
#include <type_traits>


// This is just a thin wrapper for std::mt19937_64 with the constructor that properly seeds its state from std::random_device;

class MersenneTwister
{
private:
    
    std::mt19937_64 random_engine;
    
public:
    
    using result_type = std::mt19937_64::result_type;
    using UInt        = std::mt19937_64::result_type;
    
    static constexpr std::size_t state_size = std::mt19937_64::state_size;
    
    MersenneTwister() noexcept
    {
        constexpr std::size_t seed_size
                = state_size * sizeof(result_type) / sizeof(std::random_device::result_type);
        
        std::array<std::random_device::result_type,seed_size> seed_array;
        
        std::generate( seed_array.begin(), seed_array.end(), std::random_device() );

        std::seed_seq seed ( seed_array.begin(), seed_array.end() );
        
        random_engine = std::mt19937_64( seed );
        
    }
    
    result_type operator()() noexcept
    {
        return random_engine();
    }
    
    
    
    static constexpr result_type min() noexcept
    {
        return std::mt19937_64::min();
    }
    
    static constexpr result_type max() noexcept
    {
        return std::mt19937_64::max();
    }
    
    friend bool operator ==(const MersenneTwister& lhs, const MersenneTwister& rhs) noexcept
    {
        return (lhs.random_engine == rhs.random_engine);
    }
    
    friend bool operator !=(const MersenneTwister& lhs, const MersenneTwister& rhs) noexcept
    {
        return (lhs.random_engine != rhs.random_engine);
    }
};

