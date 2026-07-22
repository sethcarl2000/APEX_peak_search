#ifndef peak_search_log_likelihood_hpp
#define peak_search_log_likelihood_hpp

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

/// @brief Computes log-likelihood for 1D histogram 
/// @param hist 1D histogram 
/// @param pdf function to return expectation value for each bin (evaluated at each bin center)
/// @return log likelihood
double log_likelihood(TH1D* hist, const std::function<double(double)>& fcn);

/// @brief Computes log-likelihood for 1D histogram 
/// @param hist vector of 1D histogram bins 
/// @param pdf function to return expectation value for each bin (evaluated at each bin center)
/// @return negative log likelihood
double negative_log_likelihood(const histo_1D_t& hist, const std::function<double(double)>& fcn);

/// @brief Computes modified NLL (combinatoric factor ignored) 
/// @param hist vector of 1D histogram bins 
/// @param pdf function to return expectation value for each bin (evaluated at each bin center)
/// @return negative log likelihood
double modified_nll(const histo_1D_t& hist, const std::function<double(double)>& fcn);

};


#endif