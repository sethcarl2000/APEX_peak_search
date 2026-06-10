#ifndef peak_search_chisquare_hpp
#define peak_search_chisquare_hpp

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

/// @brief computes chi^2 agreement between histogram and fcn  
/// @param hist 1D histogram 
/// @param pdf function to return expectation value for each bin (integrated over each bin)
/// @return chi^2
double chisquare(TH1D* hist, const std::function<double(double)>& fcn);

/// @brief computes p-value for chi-square. 
/// @param chi2 chi^2 value returned by above fcn 
/// @param n_bins number of bins in histogram
/// @return p-value (prob. that chi-square as big or larger than given val)
double chisquare_p(double chi2, int n_bins);

};


#endif