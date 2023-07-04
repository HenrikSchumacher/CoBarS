#pragma once

// SplitMix64
// Output: 64 bits
// Period: 2^64
// Footprint: 8 bytes
// Original implementation: http://prng.di.unimi.it/splitmix64.c
class SplitMix64
{
public:
    
    using UInt         = std::uint64_t;
    using state_type   = UInt;
    using result_type  = UInt;
    
    SplitMix64() noexcept
    {
        UInt32 seed [2];
        
        std::random_device rand_dev;
        
        seed[0] = rand_dev();
        seed[1] = rand_dev();
        
        state = *reinterpret_cast<UInt*>(&seed[0]);
    }
    
    explicit SplitMix64(const state_type state_ ) noexcept
    :   state(state_)
    {}
    
    result_type operator()() noexcept
    {
        UInt z = (state += static_cast<UInt>(0x9e3779b97f4a7c15));
        z = (z ^ (z >> 30)) * static_cast<UInt>(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)) * static_cast<UInt>(0x94d049bb133111eb);
        return z ^ (z >> 31);
    }
    
    template <std::size_t N>
    std::array<UInt,N> generateSeedSequence() noexcept
    {
        std::array<UInt, N> seeds = {};
        
        for( auto& seed : seeds )
        {
            seed = operator()();
        }
        
        return seeds;
    }
    
    constexpr result_type min() noexcept
    {
        return std::numeric_limits<result_type>::lowest();
    }
    
    constexpr result_type max() noexcept
    {
        return std::numeric_limits<result_type>::max();
    }
    
    state_type serialize() const noexcept
    {
        return state;
    }
    
    void deserialize(const state_type state_) noexcept
    {
        state = state_;
    }
    
    friend bool operator ==(const SplitMix64& lhs, const SplitMix64& rhs) noexcept
    {
        return (lhs.state == rhs.state);
    }
    
    friend bool operator !=(const SplitMix64& lhs, const SplitMix64& rhs) noexcept
    {
        return (lhs.state != rhs.state);
    }
    
private:
    
    state_type state;
};

