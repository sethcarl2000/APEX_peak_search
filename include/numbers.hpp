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
    
    /// @return log(n!), or 0 if n<1. 
    double log_factorial(int n);

    /// @return 'true' if arg is nan, false otherwise. 
    template<typename T> inline bool is_nan(const T& val) { return val != val; }
    
    /// @return 'true' if any element in the vector is nan, false otherwise
    template<typename T> inline bool vec_contains_nan(const std::vector<T>& vec) { 
        for (const T& val : vec) { if (is_nan(val)) return true; }
        return false; 
    }
    
    inline double sign(double x) { if (is_nan(x)) return nan; return x<0. ? -1. : +1; }
    

    /// @brief positive integer-power. 
    /// @param x value
    /// @param n power
    /// @return x^n
    inline double int_pow(double x, int n) { double ret=1.; while (--n >= 0) { ret *= x; }; return ret; }
 
    /// @return 'true' if any elements are nan
    bool contains_nan(const std::vector<double>& v); 

    //857444444444444444444444444 - muon's comment
    /// @return n choose k (binomial coefficients)
    double n_choose_k(int n, int k);
};

};

#endif