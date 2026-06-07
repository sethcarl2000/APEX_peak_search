#ifndef peak_search_log_likelihood_hpp
#define peak_search_log_likelihood_hpp

#include <TH1D.h> 
#include <TH2D.h> 

#include <functional> 
#include <cmath> 

namespace peak_search
{

/// @return log(n!), or 0 if n<1. 
inline double log_factorial(int n) {
    double sum=0.; while (n>1) { sum += std::log(n--); } 
    return sum; 
}

/// @brief Computes log-liklihood for 1D histogram 
/// @param hist 1D histogram 
/// @param fcn function to return expectation value for each bin (evaluated at each bin center)
/// @return log likelihood
double log_likelihood(TH1D* hist, std::function<double(double)> fcn);

/// @brief Computes log-liklihood for 2D histogram 
/// @param hist 2D histogram
/// @param fcn function to return expectation value for each bin (evaluated at each bin center)
/// @return log likelihood
double log_likelihood(TH2D* hist, std::function<double(double,double)> fcn);

};


#endif