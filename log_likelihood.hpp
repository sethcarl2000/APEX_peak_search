#ifndef peak_search_log_likelihood_hpp
#define peak_search_log_likelihood_hpp

#include "pdf_fcn.hpp"
#include "bininfo.hpp"

//ROOT headers
#include <TH1D.h> 
#include <TH2D.h> 

//stdlib headers 
#include <functional> 
#include <cmath> 
#include <vector> 

namespace peak_search
{

/// @return log(n!), or 0 if n<1. 
inline double log_factorial(int n) {
    double sum=0.; while (n>1) { sum += std::log(n--); } 
    return sum; 
}

/// @brief Computes log-likelihood for 1D histogram 
/// @param hist 1D histogram 
/// @param pdf function to return expectation value for each bin (evaluated at each bin center)
/// @return log likelihood
double log_likelihood(TH1D* hist, double (*fcn)(double));

/// @brief Computes log-likelihood for 1D histogram 
/// @param hist vector of 1D histogram bins 
/// @param pdf function to return expectation value for each bin (evaluated at each bin center)
/// @return negative log likelihood
double negative_log_likelihood(const histo_1D_t& hist, const std::function<double(double)>& fcn);

/// @brief Computes log-liklihood for 2D histogram 
/// @param hist 2D histogram
/// @param fcn function to return expectation value for each bin (evaluated at each bin center)
/// @return log likelihood
double log_likelihood(TH2D* hist, const PdfFcn_2D& pdf);

};


#endif