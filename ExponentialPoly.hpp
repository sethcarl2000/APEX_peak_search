#ifndef peak_search_ExponentialPoly_hpp
#define peak_search_ExponentialPoly_hpp

#include "numbers.hpp"

#include <vector> 
#include <cmath> 
#include <limits> 

namespace peak_search 
{

/// @brief Reprensents an exponentiated polynomial. uses homer technique to compute poly 
struct ExponentialPoly {
    //coefficients are stored in the following order: 
    // f(x) = std::exp( coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + coeffs[3]*x*x*x + ... )
    std::vector<double> coeffs; 

    inline double operator()(double x) const {
        if (coeffs.empty()) return numbers::nan; 
        
        double poly = coeffs.back(); 
        for (int i=((int)coeffs.size())-2; i>=0; i--) { poly = coeffs[i] + poly*x; }
        return std::exp( poly );
    }

    //coefficients are stored in the following order: 
    // f(x) = std::exp( coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + coeffs[3]*x*x*x + ... )
    inline static double Eval(double x, double *coeffs) {
        int n; 
        //check if the array of coefficients is null, or if its size is less than one (if so, return nan)
        if (!coeffs || (n = ((int)(sizeof(coeffs)/sizeof(double)))) < 1) return numbers::nan; 

        double poly = coeffs[n-1];
        for (int i=n-2; i>=0; i--) { poly = coeffs[i] + poly*x; }
        return std::exp( poly );
    }
};

}; 

#endif