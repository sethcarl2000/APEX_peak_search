#ifndef peak_search_ExponentialPoly_hpp
#define peak_search_ExponentialPoly_hpp

#include "numbers.hpp"
#include "pdf_fcn.hpp"

#include <vector> 
#include <cmath> 
#include <limits> 

namespace peak_search 
{

/// @brief Reprensents an exponentiated polynomial. uses homer technique to compute poly 
struct ExponentialPoly {
    //coefficients are stored in the following order: 
    // f(x) = std::exp( coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + coeffs[3]*x*x*x + ... )
    std::vector<double> coeffs{}; 

    double operator()(double x) const;

    //coefficients are stored in the following order: 
    // f(x) = std::exp( coeffs[0] + coeffs[1]*x + coeffs[2]*x*x + coeffs[3]*x*x*x + ... )
    inline static double Eval(double x, double *coeffs, int n);

    explicit operator PdfFcn_1D() const; 
};

static_assert(std::is_default_constructible<ExponentialPoly>::value, "ExponentialPoly struct is not default constructible"); 

}; 

#endif