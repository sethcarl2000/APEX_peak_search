#ifndef peak_search_best_likelihood_fit_hpp
#define peak_search_best_likelihood_fit_hpp

#include "FitResult.hpp"

#include <vector>
#include <TH1D.h> 

namespace peak_search
{
    
/// @brief Using the Minuit fit package wrapped with ROOT, find the best likelihood fit for the given histogram.
/// @param hist histogram to fit
/// @param fcn function to fit. subroutine will seek to match the integral of this function over each bin to the bin's contents.
/// @param params array of parameters to optimize
/// @return negative log-likelihood of best fit found. 'nan' if fit is unsuccessful. 
double best_likelihood_fit(TH1D* hist, double (*fcn)(double,const double*), std::vector<double>& params);

};
#endif