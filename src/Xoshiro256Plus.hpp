#pragma once
#include <cstdint>
#include <array>

#include "../deps/Xoshiro-cpp/XoshiroCpp.hpp"


// This is just a thin wrapper for XoshiroCpp::Xoshiro256Plus with a constructor that properly seeds its state from std::random_device;

namespace CoBarS
{
    class Xoshiro256Plus
    {
    private:
        
        XoshiroCpp::Xoshiro256Plus random_engine;
        
    public:
        
        using result_type = XoshiroCpp::Xoshiro256Plus::result_type;
        
        
        Xoshiro256Plus() noexcept
        {
            std::array<std::uint64_t, 4> seeds;
            
            std::random_device r;
            
            mptr<std::uint32_t> seeds_ = reinterpret_cast<std::uint32_t *>(&seeds[0]);
            
            for( int i = 0; i < 8; ++i )
            {
                seeds_[i] = r();
            }
            
            random_engine = XoshiroCpp::Xoshiro256Plus( seeds );
        }
        
        result_type operator()() noexcept
        {
            return random_engine();
        }
        
        
        static constexpr result_type min() noexcept
        {
            return XoshiroCpp::Xoshiro256Plus::min();
        }
        
        static constexpr result_type max() noexcept
        {
            return XoshiroCpp::Xoshiro256Plus::max();
        }
        
        std::string ClassName()
        {
            return std::string("Xoshiro256Plus");
        }
    };
    
} // namespace CoBarS
