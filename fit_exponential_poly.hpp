#ifndef peak_search_fit_exponential_poly_hpp
#define peak_search_fit_exponential_poly_hpp

#include "ExponentialPoly.hpp"
#include "bininfo.hpp"
#include "FitResult.hpp"

//ROOT headers
#include <TH1D.h> 


namespace peak_search
{

/// @brief Fit an exponentiated-polynomial to the given histogram 
/// @param hist histogram to fit 
/// @param degree degree of exponentiated polynomial 
/// @return Exponential polynomial fit-result 
FitResult<ExponentialPoly> fit_exponential_poly(TH1D* hist, int degree);

/// @brief Fit an exponentiated-polynomial to the given histogram 
/// @param hist histogram to fit 
/// @param degree degree of exponentiated polynomial 
/// @return Exponential polynomial fit-result 
FitResult<ExponentialPoly> fit_exponential_poly(histo_1D_t data, int degree);

};

#endif