#ifndef peak_search_legendre_polynomial_hpp
#define peak_search_legendre_polynomial_hpp

namespace peak_search {

/// @brief Compute legendre polynomial of order 'n' (positive integer)
/// @param x arg (interval: [0,1])
/// @param n order
/// @return legendre polynomial, normalized in standard way
double legendre_polynomial(double x, int n); 

/// @brief First deriv of legendre polynomial: (d/dx) P_n(x) 
double legendre_polynomial_dx(double x, int n); 

/// @brief Second deriv of legendre polynomial: (d^2/dx^2) P_n(x) 
double legendre_polynomial_ddx(double x, int n); 

};

#endif