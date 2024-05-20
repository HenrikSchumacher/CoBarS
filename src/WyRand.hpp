#pragma once
#include <cstdint>
#include <array>

#include "../submodules/wy/wy.hpp"


// This is just a thin wrapper for wy::rand with a constructor that properly seeds its state from std::random_device;

namespace CoBarS
{
    
    class WyRand
    {
    public:
        
        using result_type = wy::rand::result_type;
    
        
    private:
        
        wy::rand random_engine;
        
    public:
        
        WyRand() 
        {}
        
        ~WyRand() = default;

        result_type operator()() noexcept
        {
            return random_engine();
        }
        
        static constexpr result_type min() noexcept
        {
            return wy::rand::min();
        }
        
        static constexpr result_type max() noexcept
        {
            return wy::rand::max();
        }
        
//        friend bool operator ==(const WyRand& lhs, const WyRand& rhs) noexcept
//        {
//            return (lhs.random_engine == rhs.random_engine);
//        }
//        
//        friend bool operator !=(const WyRand& lhs, const WyRand& rhs) noexcept
//        {
//            return (lhs.random_engine != rhs.random_engine);
//        }
        
        std::string ClassName()
        {
            return std::string("WyRand");
        }
    };
    
} // namespace CoBarS

