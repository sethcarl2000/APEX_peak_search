#ifndef peak_search_numbers_hpp
#define peak_search_numbers_hpp

#include <limits>
#include <vector> 

namespace peak_search
{
//be careful; this will clash with the std::numbers namespace if you are working in both namespaces. 
namespace numbers 
{
    constexpr double pi = 3.1415926536;

    constexpr double nan = std::numeric_limits<double>::quiet_NaN(); 

    /// @return 'true' if arg is nan, false otherwise. 
    inline bool is_nan(double x) { return x != x; }

    /// @brief positive integer-power. 
    /// @param x value
    /// @param n power
    /// @return x^n
    inline double int_pow(double x, int n) { double ret=x; while (--n > 0) { ret *= x; }; return ret; }
 
    /// @return 'true' if any elements are nan
    bool contains_nan(const std::vector<double>& v); 
};

};

#endif